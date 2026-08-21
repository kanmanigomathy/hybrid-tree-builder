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


# ============================================================
# SEQUENCE TYPE
# ============================================================

def detect_sequence_type(sequence):

    sequence = sequence.upper().replace("-", "")
    chars = set(sequence)

    dna_chars = set("ACGTN")
    rna_chars = set("ACGUN")

    if chars and chars.issubset(dna_chars):
        return "DNA"

    if chars and chars.issubset(rna_chars):
        return "RNA"

    return "Protein"


# ============================================================
# GET UPLOADED FILES
# ============================================================

def get_uploaded_files():

    files = request.files.getlist("files")

    # Compatibility with old single-file upload
    if not files:
        files = request.files.getlist("file")

    return [
        f for f in files
        if f and f.filename
    ]


# ============================================================
# PARSE MULTIPLE FASTA FILES
# ============================================================

def parse_fasta_files(files):

    sequences = []
    errors = []

    if not files:
        return [], ["Please upload at least one FASTA file."]

    for file in files:

        filename = file.filename

        if not filename.lower().endswith(
            (".fasta", ".fa", ".fas")
        ):
            errors.append(
                f"{filename}: unsupported file type."
            )
            continue

        try:

            raw = file.read()

            text = raw.decode("utf-8-sig")

            records = list(
                SeqIO.parse(
                    StringIO(text),
                    "fasta"
                )
            )

            if not records:

                errors.append(
                    f"{filename}: no FASTA sequences found."
                )

                continue

            for record in records:

                sequence = (
                    str(record.seq)
                    .upper()
                    .replace(" ", "")
                    .replace("\r", "")
                    .replace("\n", "")
                )

                if not sequence:

                    errors.append(
                        f"{filename} / {record.id}: empty sequence."
                    )

                    continue

                sequences.append({

                    "name": record.id,

                    "length": len(sequence),

                    "sequence": sequence,

                    "type": detect_sequence_type(
                        sequence
                    ),

                    "source_file": filename

                })

        except UnicodeDecodeError:

            errors.append(
                f"{filename}: invalid encoding. "
                f"Save the FASTA file as UTF-8."
            )

        except Exception as e:

            errors.append(
                f"{filename}: invalid FASTA format: {e}"
            )

    return sequences, errors


# ============================================================
# VALIDATION STATISTICS
# ============================================================

def validation_statistics(sequences):

    lengths = [
        s["length"]
        for s in sequences
    ]

    types = [
        s["type"]
        for s in sequences
    ]

    if types and all(
        t == types[0]
        for t in types
    ):
        overall_type = types[0]
    else:
        overall_type = "Mixed"

    return {

        "number_of_sequences":
            len(sequences),

        "sequence_type":
            overall_type,

        "minimum_length":
            min(lengths)
            if lengths else 0,

        "maximum_length":
            max(lengths)
            if lengths else 0,

        "average_length":
            round(
                statistics.mean(lengths),
                2
            )
            if lengths else 0,

        "total_length":
            sum(lengths),

        "number_of_files":
            len(
                set(
                    s["source_file"]
                    for s in sequences
                )
            )
    }


# ============================================================
# HOME
# ============================================================

@app.route("/")
def home():

    return render_template(
        "index.html"
    )


# ============================================================
# VALIDATE
# ============================================================

@app.route(
    "/validate",
    methods=["POST"]
)
def validate():

    files = get_uploaded_files()

    if not files:

        return jsonify({

            "success": False,

            "message":
                "Please upload at least one FASTA file."

        }), 400

    sequences, errors = parse_fasta_files(
        files
    )

    if not sequences:

        return jsonify({

            "success": False,

            "errors":
                errors

        }), 400

    return jsonify({

        "success": True,

        "message":
            "FASTA validation successful.",

        "statistics":
            validation_statistics(
                sequences
            ),

        "sequences":
            sequences,

        "warnings":
            errors

    })


# ============================================================
# PAIRWISE ALIGNMENT
# ============================================================

def pairwise_to_gapped(
    reference,
    sequence
):

    aligner = PairwiseAligner()

    aligner.match_score = 2

    aligner.mismatch_score = -1

    aligner.open_gap_score = -2

    aligner.extend_gap_score = -0.5

    alignment = aligner.align(
        reference,
        sequence
    )[0]

    ref_output = []

    seq_output = []

    coordinates = alignment.coordinates

    for i in range(
        coordinates.shape[1] - 1
    ):

        r0 = int(
            coordinates[0, i]
        )

        r1 = int(
            coordinates[0, i + 1]
        )

        s0 = int(
            coordinates[1, i]
        )

        s1 = int(
            coordinates[1, i + 1]
        )

        dr = r1 - r0
        ds = s1 - s0

        if dr and ds:

            ref_output.append(
                reference[r0:r1]
            )

            seq_output.append(
                sequence[s0:s1]
            )

        elif dr and not ds:

            ref_output.append(
                reference[r0:r1]
            )

            seq_output.append(
                "-" * dr
            )

        elif ds and not dr:

            ref_output.append(
                "-" * ds
            )

            seq_output.append(
                sequence[s0:s1]
            )

    return (
        "".join(ref_output),
        "".join(seq_output)
    )


# ============================================================
# PROGRESSIVE ALIGNMENT
# ============================================================

def merge_alignment(
    rows,
    new_reference,
    new_sequence
):

    old_reference = rows[0]["sequence"]

    if old_reference == new_reference:

        return rows + [{

            "name": "",
            "sequence": new_sequence

        }]

    merged = [

        {
            "name": row["name"],
            "sequence": ""
        }

        for row in rows

    ]

    new_row = ""

    i = 0
    j = 0

    while (
        i < len(old_reference)
        or
        j < len(new_reference)
    ):

        old_char = (
            old_reference[i]
            if i < len(old_reference)
            else None
        )

        new_char = (
            new_reference[j]
            if j < len(new_reference)
            else None
        )

        if (
            old_char == new_char
        ):

            for k, row in enumerate(rows):

                merged[k]["sequence"] += (
                    row["sequence"][i]
                )

            new_row += new_sequence[j]

            i += 1
            j += 1

        elif (
            old_char == "-"
            and
            new_char != "-"
        ):

            for k, row in enumerate(rows):

                merged[k]["sequence"] += (
                    row["sequence"][i]
                )

            new_row += "-"

            i += 1

        elif (
            new_char == "-"
            and
            old_char != "-"
        ):

            for k in range(len(rows)):

                merged[k]["sequence"] += "-"

            new_row += new_sequence[j]

            j += 1

        else:

            if old_char is not None:

                for k, row in enumerate(rows):

                    merged[k]["sequence"] += (
                        row["sequence"][i]
                    )

                new_row += "-"

                i += 1

            if new_char is not None:

                for k in range(len(rows)):

                    merged[k]["sequence"] += "-"

                new_row += new_sequence[j]

                j += 1

    merged.append({

        "name": "",

        "sequence": new_row

    })

    return merged


# ============================================================
# CREATE MULTIPLE SEQUENCE ALIGNMENT
# ============================================================

def create_alignment(sequences):

    if len(sequences) < 2:

        return (
            None,
            "At least two sequences are required."
        )

    rows = [{

        "name":
            sequences[0]["name"],

        "sequence":
            sequences[0]["sequence"]

    }]

    for item in sequences[1:]:

        reference = (
            rows[0]["sequence"]
            .replace("-", "")
        )

        new_reference, new_sequence = (
            pairwise_to_gapped(
                reference,
                item["sequence"]
            )
        )

        rows = merge_alignment(
            rows,
            new_reference,
            new_sequence
        )

        rows[-1]["name"] = item["name"]

    # Restore names
    for index, row in enumerate(rows):

        row["name"] = sequences[index]["name"]

    # Same alignment length
    max_length = max(
        len(row["sequence"])
        for row in rows
    )

    for row in rows:

        row["sequence"] = (
            row["sequence"]
            .ljust(
                max_length,
                "-"
            )
        )

    return rows, None


# ============================================================
# ALIGN API
# ============================================================

@app.route(
    "/align",
    methods=["POST"]
)
def align_sequences():

    start_time = time.time()

    files = get_uploaded_files()

    if not files:

        return jsonify({

            "success": False,

            "message":
                "Please upload at least one FASTA file."

        }), 400

    sequences, errors = parse_fasta_files(
        files
    )

    if not sequences:

        return jsonify({

            "success": False,

            "message":
                "No valid sequences found.",

            "errors":
                errors

        }), 400

    aligned, error = create_alignment(
        sequences
    )

    if error:

        return jsonify({

            "success": False,

            "message":
                error

        }), 400

    execution_time = round(
        time.time() - start_time,
        4
    )

    return jsonify({

        "success": True,

        "message":
            "Multiple sequence alignment completed.",

        "number_of_sequences":
            len(aligned),

        "alignment_length":
            len(aligned[0]["sequence"]),

        "execution_time":
            execution_time,

        "alignment":
            aligned,

        "warnings":
            errors

    })


# ============================================================
# NEIGHBOR JOINING TREE
# ============================================================

def build_neighbor_joining_tree(
    aligned_sequences
):

    if len(aligned_sequences) < 3:

        return (
            None,
            "At least 3 sequences are required."
        )

    records = [

        SeqRecord(
            Seq(item["sequence"]),
            id=item["name"],
            description=""
        )

        for item in aligned_sequences

    ]

    alignment = MultipleSeqAlignment(
        records
    )

    calculator = DistanceCalculator(
        "identity"
    )

    distance_matrix = (
        calculator.get_distance(
            alignment
        )
    )

    constructor = (
        DistanceTreeConstructor()
    )

    tree = constructor.nj(
        distance_matrix
    )

    try:

        tree.root_at_midpoint()

    except Exception:

        pass

    output = StringIO()

    Phylo.write(
        tree,
        output,
        "newick"
    )

    return (
        output.getvalue(),
        None
    )


# ============================================================
# TREE API
# ============================================================

@app.route(
    "/build-tree",
    methods=["POST"]
)
def build_tree():

    start_time = time.time()

    files = get_uploaded_files()

    if not files:

        return jsonify({

            "success": False,

            "message":
                "Please upload FASTA files."

        }), 400

    sequences, errors = parse_fasta_files(
        files
    )

    if not sequences:

        return jsonify({

            "success": False,

            "errors":
                errors

        }), 400

    aligned, alignment_error = (
        create_alignment(
            sequences
        )
    )

    if alignment_error:

        return jsonify({

            "success": False,

            "message":
                alignment_error

        }), 400

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

        }), 400

    return jsonify({

        "success": True,

        "method":
            "Neighbor Joining",

        "number_of_sequences":
            len(aligned),

        "execution_time":
            round(
                time.time() - start_time,
                4
            ),

        "newick":
            newick,

        "warnings":
            errors

    })


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    app.run(
        host="127.0.0.1",
        port=5000,
        debug=True
    )