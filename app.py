from flask import Flask, render_template, request, jsonify
from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from io import StringIO
import statistics
import time


app = Flask(__name__)


# =========================================================
# SEQUENCE TYPE DETECTION
# =========================================================

def detect_sequence_type(sequence):

    sequence = sequence.upper()

    dna_chars = set("ACGTN")
    rna_chars = set("ACGUN")

    sequence_chars = set(sequence)

    if sequence_chars.issubset(dna_chars):
        return "DNA"

    if sequence_chars.issubset(rna_chars):
        return "RNA"

    return "Protein"


# =========================================================
# FASTA VALIDATION
# =========================================================

def validate_fasta(text):

    sequences = []
    errors = []

    try:

        records = list(
            SeqIO.parse(
                StringIO(text),
                "fasta"
            )
        )

        if len(records) == 0:

            errors.append(
                "No FASTA sequences found."
            )

            return sequences, errors


        for record in records:

            sequence = str(
                record.seq
            ).upper()


            if len(sequence) == 0:

                errors.append(
                    f"{record.id}: Empty sequence."
                )

                continue


            sequence_type = detect_sequence_type(
                sequence
            )


            sequence_data = {

                "name": record.id,

                "length": len(sequence),

                "sequence": sequence,

                "type": sequence_type
            }


            sequences.append(
                sequence_data
            )


    except Exception as error:

        errors.append(
            "Invalid FASTA format: "
            + str(error)
        )


    return sequences, errors


# =========================================================
# HOME PAGE
# =========================================================

@app.route("/")
def home():

    return render_template(
        "index.html"
    )


# =========================================================
# FASTA VALIDATION API
# =========================================================

@app.route(
    "/validate",
    methods=["POST"]
)
def validate():

    if "file" not in request.files:

        return jsonify({

            "success": False,

            "message":
                "Please upload a FASTA file."
        })


    file = request.files["file"]


    if file.filename == "":

        return jsonify({

            "success": False,

            "message":
                "No file selected."
        })


    try:

        text = file.read().decode(
            "utf-8"
        )

    except UnicodeDecodeError:

        return jsonify({

            "success": False,

            "message":
                "Invalid text encoding."
        })


    sequences, errors = validate_fasta(
        text
    )


    if len(errors) > 0:

        return jsonify({

            "success": False,

            "errors": errors
        })


    if len(sequences) == 0:

        return jsonify({

            "success": False,

            "message":
                "No valid sequences found."
        })


    lengths = [

        item["length"]

        for item in sequences
    ]


    sequence_types = [

        item["type"]

        for item in sequences
    ]


    if all(
        item == sequence_types[0]
        for item in sequence_types
    ):

        overall_type = sequence_types[0]

    else:

        overall_type = "Mixed"


    statistics_data = {

        "number_of_sequences":
            len(sequences),

        "sequence_type":
            overall_type,

        "minimum_length":
            min(lengths),

        "maximum_length":
            max(lengths),

        "average_length":
            round(
                statistics.mean(lengths),
                2
            )
    }


    return jsonify({

        "success": True,

        "message":
            "FASTA validation successful.",

        "statistics":
            statistics_data,

        "sequences":
            sequences
    })


# =========================================================
# MULTIPLE SEQUENCE ALIGNMENT
# =========================================================

def create_alignment(sequences):

    if len(sequences) < 2:

        return (
            None,
            "At least two sequences are required."
        )


    reference = sequences[0]["sequence"]


    aligner = PairwiseAligner()

    aligner.match_score = 2

    aligner.mismatch_score = -1

    aligner.open_gap_score = -2

    aligner.extend_gap_score = -0.5


    aligned_sequences = []


    aligned_sequences.append({

        "name":
            sequences[0]["name"],

        "sequence":
            reference
    })


    for item in sequences[1:]:

        sequence = item["sequence"]


        alignment = aligner.align(
            reference,
            sequence
        )[0]


        aligned_reference = []

        aligned_sequence = []


        ref_blocks = alignment.aligned[0]

        seq_blocks = alignment.aligned[1]


        ref_position = 0

        seq_position = 0


        for i in range(
            len(ref_blocks)
        ):

            ref_start = ref_blocks[i][0]

            ref_end = ref_blocks[i][1]

            seq_start = seq_blocks[i][0]

            seq_end = seq_blocks[i][1]


            # Reference gap

            if ref_start > ref_position:

                gap_length = (
                    ref_start -
                    ref_position
                )


                aligned_reference.extend(

                    reference[
                        ref_position:
                        ref_start
                    ]
                )


                aligned_sequence.extend(

                    "-" * gap_length
                )


            # Sequence gap

            if seq_start > seq_position:

                gap_length = (
                    seq_start -
                    seq_position
                )


                aligned_reference.extend(

                    "-" * gap_length
                )


                aligned_sequence.extend(

                    sequence[
                        seq_position:
                        seq_start
                    ]
                )


            # Matching block

            aligned_reference.extend(

                reference[
                    ref_start:
                    ref_end
                ]
            )


            aligned_sequence.extend(

                sequence[
                    seq_start:
                    seq_end
                ]
            )


            ref_position = ref_end

            seq_position = seq_end


        # Remaining reference

        if ref_position < len(reference):

            aligned_reference.extend(

                reference[
                    ref_position:
                ]
            )


            aligned_sequence.extend(

                "-" *
                (
                    len(reference)
                    - ref_position
                )
            )


        # Remaining sequence

        if seq_position < len(sequence):

            aligned_reference.extend(

                "-" *
                (
                    len(sequence)
                    - seq_position
                )
            )


            aligned_sequence.extend(

                sequence[
                    seq_position:
                ]
            )


        aligned_sequences.append({

            "name":
                item["name"],

            "sequence":
                "".join(
                    aligned_sequence
                )
        })


    # Make all sequences equal length

    max_length = max(

        len(item["sequence"])

        for item in aligned_sequences
    )


    for item in aligned_sequences:

        current_length = len(
            item["sequence"]
        )


        if current_length < max_length:

            item["sequence"] += (

                "-" *
                (
                    max_length -
                    current_length
                )
            )


    return aligned_sequences, None


# =========================================================
# ALIGNMENT API
# =========================================================

@app.route(
    "/align",
    methods=["POST"]
)
def align_sequences():

    start_time = time.time()


    if "file" not in request.files:

        return jsonify({

            "success": False,

            "message":
                "Please upload a FASTA file."
        })


    file = request.files["file"]


    try:

        text = file.read().decode(
            "utf-8"
        )

    except UnicodeDecodeError:

        return jsonify({

            "success": False,

            "message":
                "Invalid FASTA encoding."
        })


    sequences, errors = validate_fasta(
        text
    )


    if len(errors) > 0:

        return jsonify({

            "success": False,

            "errors": errors
        })


    aligned, alignment_error = create_alignment(
        sequences
    )


    if alignment_error:

        return jsonify({

            "success": False,

            "message":
                alignment_error
        })


    alignment_length = len(

        aligned[0]["sequence"]
    )


    execution_time = round(

        time.time() -
        start_time,

        4
    )


    return jsonify({

        "success": True,

        "message":
            "Multiple sequence alignment completed.",

        "number_of_sequences":
            len(aligned),

        "alignment_length":
            alignment_length,

        "execution_time":
            execution_time,

        "alignment":
            aligned
    })


# =========================================================
# NEIGHBOR JOINING TREE
# =========================================================

def build_neighbor_joining_tree(
    aligned_sequences
):

    if len(aligned_sequences) < 3:

        return (
            None,
            "At least 3 sequences are required to build a phylogenetic tree."
        )


    records = []


    for item in aligned_sequences:

        records.append(

            SeqRecord(

                Seq(
                    item["sequence"]
                ),

                id=item["name"],

                description=""
            )
        )


    alignment = MultipleSeqAlignment(
        records
    )


    # Calculate evolutionary distance

    calculator = DistanceCalculator(
        "identity"
    )


    distance_matrix = calculator.get_distance(
        alignment
    )


    # Neighbor Joining

    constructor = DistanceTreeConstructor()


    tree = constructor.nj(
        distance_matrix
    )


    # Midpoint rooting

    try:

        tree.root_at_midpoint()

    except Exception:

        pass


    # Convert tree to Newick

    output = StringIO()


    Phylo.write(

        tree,

        output,

        "newick"
    )


    newick = output.getvalue()


    return newick, None


# =========================================================
# NEIGHBOR JOINING API
# =========================================================

@app.route(
    "/build-tree",
    methods=["POST"]
)
def build_tree():

    start_time = time.time()


    if "file" not in request.files:

        return jsonify({

            "success": False,

            "message":
                "Please upload a FASTA file."
        })


    file = request.files["file"]


    try:

        text = file.read().decode(
            "utf-8"
        )

    except UnicodeDecodeError:

        return jsonify({

            "success": False,

            "message":
                "Invalid FASTA encoding."
        })


    # Validate

    sequences, errors = validate_fasta(
        text
    )


    if len(errors) > 0:

        return jsonify({

            "success": False,

            "errors": errors
        })


    # Alignment

    aligned, alignment_error = create_alignment(
        sequences
    )


    if alignment_error:

        return jsonify({

            "success": False,

            "message":
                alignment_error
        })


    # NJ tree

    newick, tree_error = (
        build_neighbor_joining_tree(
            aligned
        )
    )


    if tree_error:

        return jsonify({

            "success": False,

            "message":
                tree_error
        })


    execution_time = round(

        time.time() -
        start_time,

        4
    )


    return jsonify({

        "success": True,

        "method":
            "Neighbor Joining",

        "message":
            "Neighbor Joining tree constructed successfully.",

        "number_of_sequences":
            len(aligned),

        "execution_time":
            execution_time,

        "newick":
            newick
    })


# =========================================================
# START SERVER
# =========================================================

if __name__ == "__main__":

    app.run(

        host="127.0.0.1",

        port=5000,

        debug=True
    )