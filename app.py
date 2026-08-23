from flask import Flask, render_template, request, jsonify
from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import Counter
from io import StringIO
import random
import statistics
import tempfile
import subprocess
import os
import time
import uuid


app = Flask(__name__)

app.config["MAX_CONTENT_LENGTH"] = 25 * 1024 * 1024

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FASTTREE_PATH = os.path.join(BASE_DIR, "FastTree.exe")

# Stores adaptive-bootstrap sessions temporarily
BOOTSTRAP_JOBS = {}


# ============================================================
# FASTA PROCESSING
# ============================================================

def detect_sequence_type(sequence):

    sequence = sequence.upper().replace("-", "")

    chars = set(sequence)

    if chars and chars.issubset(set("ACGTN")):
        return "DNA"

    if chars and chars.issubset(set("ACGUN")):
        return "RNA"

    return "Protein"


def get_uploaded_files():

    files = request.files.getlist("files")

    if not files:
        files = request.files.getlist("file")

    return [
        f for f in files
        if f and f.filename
    ]


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

            text = file.read().decode("utf-8-sig")

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
                "Save FASTA as UTF-8."
            )

        except Exception as exc:

            errors.append(
                f"{filename}: invalid FASTA format: {exc}"
            )

    return sequences, errors


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
            min(lengths) if lengths else 0,

        "maximum_length":
            max(lengths) if lengths else 0,

        "average_length":
            round(
                statistics.mean(lengths),
                2
            ) if lengths else 0,

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
# MULTIPLE SEQUENCE ALIGNMENT
# ============================================================

def pairwise_to_gapped(reference, sequence):

    aligner = PairwiseAligner()

    aligner.match_score = 2

    aligner.mismatch_score = -1

    aligner.open_gap_score = -2

    aligner.extend_gap_score = -0.5

    alignment = aligner.align(
        reference,
        sequence
    )[0]

    coordinates = alignment.coordinates

    ref_output = []

    seq_output = []

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


def merge_alignment(
    rows,
    new_reference,
    new_sequence
):

    old_reference = rows[0]["sequence"]

    if old_reference == new_reference:

        return rows + [
            {
                "name": "",
                "sequence": new_sequence
            }
        ]

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

        if old_char == new_char:

            for k, row in enumerate(rows):

                merged[k]["sequence"] += (
                    row["sequence"][i]
                )

            new_row += new_sequence[j]

            i += 1
            j += 1

        elif (
            old_char == "-"
            and new_char != "-"
        ):

            for k, row in enumerate(rows):

                merged[k]["sequence"] += (
                    row["sequence"][i]
                )

            new_row += "-"

            i += 1

        elif (
            new_char == "-"
            and old_char != "-"
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


def create_alignment(sequences):

    if len(sequences) < 2:

        return (
            None,
            "At least two sequences are required."
        )

    rows = [

        {
            "name":
                sequences[0]["name"],

            "sequence":
                sequences[0]["sequence"]
        }

    ]

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

    for index, row in enumerate(rows):

        row["name"] = sequences[index]["name"]

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


def write_alignment_file(
    aligned_sequences,
    file_path
):

    records = [

        SeqRecord(
            Seq(item["sequence"]),
            id=item["name"],
            description=""
        )

        for item in aligned_sequences
    ]

    with open(
        file_path,
        "w",
        encoding="utf-8"
    ) as handle:

        SeqIO.write(
            MultipleSeqAlignment(records),
            handle,
            "fasta"
        )


# ============================================================
# TREE CONSTRUCTION
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

    return output.getvalue(), None


def build_ml_tree(
    alignment_file,
    sequence_type,
    starting_tree=None
):

    if not os.path.exists(
        FASTTREE_PATH
    ):

        return (
            None,
            "FastTree.exe was not found. "
            "Place FastTree.exe beside app.py."
        )

    command = [
        FASTTREE_PATH
    ]

    if sequence_type in (
        "DNA",
        "RNA"
    ):

        command.append("-nt")

    temp_tree_file = None

    try:

        if starting_tree:

            handle = tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".nwk",
                delete=False,
                encoding="utf-8"
            )

            handle.write(
                starting_tree
            )

            handle.close()

            temp_tree_file = handle.name

            command.extend([
                "-intree",
                temp_tree_file
            ])

        command.append(
            alignment_file
        )

        result = subprocess.run(

            command,

            capture_output=True,

            text=True,

            timeout=180
        )

        if result.returncode != 0:

            return (
                None,
                result.stderr.strip()
                or "FastTree failed."
            )

        newick = (
            result.stdout.strip()
        )

        if not newick:

            return (
                None,
                "FastTree returned an empty tree."
            )

        return newick, None

    except subprocess.TimeoutExpired:

        return (
            None,
            "FastTree took too long."
        )

    except Exception as exc:

        return None, str(exc)

    finally:

        if (
            temp_tree_file
            and
            os.path.exists(
                temp_tree_file
            )
        ):

            try:
                os.remove(
                    temp_tree_file
                )
            except OSError:
                pass


def prepare_alignment_and_type(
    files
):

    sequences, errors = (
        parse_fasta_files(files)
    )

    if not sequences:

        return (
            None,
            None,
            None,
            errors,
            "No valid sequences found."
        )

    aligned, alignment_error = (
        create_alignment(
            sequences
        )
    )

    if alignment_error:

        return (
            None,
            None,
            None,
            errors,
            alignment_error
        )

    if len(aligned) < 3:

        return (
            None,
            None,
            None,
            errors,
            "At least 3 sequences are required."
        )

    types = [
        s["type"]
        for s in sequences
    ]

    if (
        types
        and
        all(
            t == types[0]
            for t in types
        )
    ):

        sequence_type = types[0]

    else:

        sequence_type = "Mixed"

    return (
        sequences,
        aligned,
        sequence_type,
        errors,
        None
    )


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

    sequences, errors = (
        parse_fasta_files(files)
    )

    if not sequences:

        return jsonify({

            "success": False,

            "message":
                "No valid sequences found.",

            "errors": errors

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
# ALIGN
# ============================================================

@app.route(
    "/align",
    methods=["POST"]
)
def align_sequences():

    start = time.time()

    files = get_uploaded_files()

    if not files:

        return jsonify({

            "success": False,

            "message":
                "Please upload FASTA files."

        }), 400

    sequences, errors = (
        parse_fasta_files(files)
    )

    if not sequences:

        return jsonify({

            "success": False,

            "message":
                "No valid sequences found.",

            "errors":
                errors

        }), 400

    aligned, error = (
        create_alignment(
            sequences
        )
    )

    if error:

        return jsonify({

            "success": False,

            "message": error

        }), 400

    return jsonify({

        "success": True,

        "message":
            "Multiple sequence alignment completed.",

        "number_of_sequences":
            len(aligned),

        "alignment_length":
            len(
                aligned[0]["sequence"]
            ),

        "execution_time":
            round(
                time.time() - start,
                4
            ),

        "alignment":
            aligned,

        "warnings":
            errors

    })


# ============================================================
# NJ
# ============================================================

@app.route(
    "/build-nj",
    methods=["POST"]
)
def build_nj():

    start = time.time()

    files = get_uploaded_files()

    (
        sequences,
        aligned,
        sequence_type,
        errors,
        error
    ) = prepare_alignment_and_type(
        files
    )

    if error:

        return jsonify({

            "success": False,

            "message": error,

            "warnings": errors

        }), 400

    nj_newick, nj_error = (
        build_neighbor_joining_tree(
            aligned
        )
    )

    if nj_error:

        return jsonify({

            "success": False,

            "message": nj_error

        }), 400

    return jsonify({

        "success": True,

        "tree_type":
            "Neighbor Joining",

        "nj_newick":
            nj_newick,

        "number_of_sequences":
            len(aligned),

        "alignment_length":
            len(
                aligned[0]["sequence"]
            ),

        "sequence_type":
            sequence_type,

        "execution_time":
            round(
                time.time() - start,
                4
            ),

        "warnings":
            errors

    })


# ============================================================
# ML
# ============================================================

@app.route(
    "/build-ml",
    methods=["POST"]
)
def build_ml():

    start = time.time()

    files = get_uploaded_files()

    (
        sequences,
        aligned,
        sequence_type,
        errors,
        error
    ) = prepare_alignment_and_type(
        files
    )

    if error:

        return jsonify({

            "success": False,

            "message": error

        }), 400

    temp_alignment = None

    try:

        handle = tempfile.NamedTemporaryFile(

            mode="w",

            suffix=".fasta",

            delete=False,

            encoding="utf-8"
        )

        temp_alignment = handle.name

        handle.close()

        write_alignment_file(
            aligned,
            temp_alignment
        )

        ml_newick, ml_error = (
            build_ml_tree(
                temp_alignment,
                sequence_type
            )
        )

        if ml_error:

            return jsonify({

                "success": False,

                "message": ml_error

            }), 500

        return jsonify({

            "success": True,

            "tree_type":
                "Maximum Likelihood",

            "ml_newick":
                ml_newick,

            "number_of_sequences":
                len(aligned),

            "alignment_length":
                len(
                    aligned[0]["sequence"]
                ),

            "sequence_type":
                sequence_type,

            "execution_time":
                round(
                    time.time() - start,
                    4
                ),

            "warnings":
                errors

        })

    finally:

        if (
            temp_alignment
            and
            os.path.exists(
                temp_alignment
            )
        ):

            try:
                os.remove(
                    temp_alignment
                )
            except OSError:
                pass


# ============================================================
# HYBRID TREE
# ============================================================

@app.route(
    "/build-hybrid",
    methods=["POST"]
)
def build_hybrid():

    start = time.time()

    files = get_uploaded_files()

    (
        sequences,
        aligned,
        sequence_type,
        errors,
        error
    ) = prepare_alignment_and_type(
        files
    )

    if error:

        return jsonify({

            "success": False,

            "message": error

        }), 400

    # STEP 1: NJ
    nj_newick, nj_error = (
        build_neighbor_joining_tree(
            aligned
        )
    )

    if nj_error:

        return jsonify({

            "success": False,

            "message": nj_error

        }), 400

    temp_alignment = None

    try:

        handle = tempfile.NamedTemporaryFile(

            mode="w",

            suffix=".fasta",

            delete=False,

            encoding="utf-8"
        )

        temp_alignment = handle.name

        handle.close()

        write_alignment_file(
            aligned,
            temp_alignment
        )

        # STEP 2: ML refinement
        ml_newick, ml_error = (
            build_ml_tree(

                temp_alignment,

                sequence_type,

                starting_tree=nj_newick
            )
        )

        if ml_error:

            return jsonify({

                "success": False,

                "message":
                    "NJ succeeded but ML refinement failed.",

                "nj_newick":
                    nj_newick,

                "ml_error":
                    ml_error

            }), 500

        return jsonify({

            "success": True,

            "method":
                "Hybrid: NJ starting topology + ML refinement",

            "nj_newick":
                nj_newick,

            "ml_newick":
                ml_newick,

            "final_newick":
                ml_newick,

            "newick":
                ml_newick,

            "number_of_sequences":
                len(aligned),

            "alignment_length":
                len(
                    aligned[0]["sequence"]
                ),

            "sequence_type":
                sequence_type,

            "execution_time":
                round(
                    time.time() - start,
                    4
                ),

            "warnings":
                errors

        })

    finally:

        if (
            temp_alignment
            and
            os.path.exists(
                temp_alignment
            )
        ):

            try:
                os.remove(
                    temp_alignment
                )
            except OSError:
                pass


# ============================================================
# BOOTSTRAP
# ============================================================

def resample_alignment(
    aligned_sequences
):

    length = len(
        aligned_sequences[0]["sequence"]
    )

    positions = [

        random.randrange(length)

        for _ in range(length)

    ]

    return [

        {

            "name":
                item["name"],

            "sequence":
                "".join(
                    item["sequence"][p]
                    for p in positions
                )

        }

        for item in aligned_sequences
    ]


def get_tree_clades(newick):

    try:

        tree = Phylo.read(
            StringIO(newick),
            "newick"
        )

        terminals = len(
            tree.get_terminals()
        )

        clades = []

        for clade in tree.get_nonterminals():

            names = tuple(
                sorted(
                    terminal.name
                    for terminal
                    in clade.get_terminals()
                    if terminal.name
                )
            )

            if (
                1 < len(names)
                < terminals
            ):

                clades.append(names)

        return set(clades)

    except Exception:

        return set()


def calculate_support(
    original_clades,
    replicate_clades
):

    total = len(
        replicate_clades
    )

    if total == 0:

        return {
            clade: 0.0
            for clade
            in original_clades
        }

    support = {}

    for clade in original_clades:

        count = sum(

            1

            for rep
            in replicate_clades

            if clade in rep
        )

        support[clade] = round(

            (
                count / total
            ) * 100,

            2

        )

    return support


def branch_results_from_support(
    support
):

    result = []

    sorted_items = sorted(

        support.items(),

        key=lambda x: x[1],

        reverse=True

    )

    for index, (
        clade,
        value
    ) in enumerate(
        sorted_items,
        start=1
    ):

        if value >= 80:

            confidence = "High"

        elif value >= 50:

            confidence = "Moderate"

        else:

            confidence = "Low"

        result.append({

            "id":
                index,

            "branch":
                " + ".join(clade),

            "members":
                list(clade),

            "support":
                value,

            "confidence":
                confidence

        })

    return result


def confidence_summary(
    support
):

    values = list(
        support.values()
    )

    if not values:

        return {

            "average_confidence":
                0,

            "high_confidence":
                0,

            "moderate_confidence":
                0,

            "low_confidence":
                0
        }

    return {

        "average_confidence":
            round(
                statistics.mean(values),
                2
            ),

        "high_confidence":
            sum(
                v >= 80
                for v in values
            ),

        "moderate_confidence":
            sum(
                50 <= v < 80
                for v in values
            ),

        "low_confidence":
            sum(
                v < 50
                for v in values
            )
    }


def run_bootstrap_batch(

    aligned,

    original_clades,

    replicate_clades,

    number

):

    for _ in range(number):

        replicate = (
            resample_alignment(
                aligned
            )
        )

        newick, error = (
            build_neighbor_joining_tree(
                replicate
            )
        )

        if not error:

            replicate_clades.append(

                get_tree_clades(
                    newick
                )
            )

    return calculate_support(

        original_clades,

        replicate_clades
    )


# ============================================================
# BOOTSTRAP ROUND 1
# ============================================================

@app.route(
    "/bootstrap-round1",
    methods=["POST"]
)
def bootstrap_round1():

    start = time.time()

    files = get_uploaded_files()

    try:

        initial = int(
            request.form.get(
                "initial_replicates",
                100
            )
        )

        batch = int(
            request.form.get(
                "batch_size",
                100
            )
        )

        maximum = int(
            request.form.get(
                "max_replicates",
                1000
            )
        )

        threshold = float(
            request.form.get(
                "stability_threshold",
                2
            )
        )

    except ValueError:

        return jsonify({

            "success": False,

            "message":
                "Invalid bootstrap parameters."

        }), 400

    initial = max(
        10,
        min(initial, 1000)
    )

    batch = max(
        10,
        min(batch, 500)
    )

    maximum = max(
        initial,
        min(maximum, 2000)
    )

    threshold = max(
        0.1,
        min(threshold, 20)
    )

    (
        sequences,
        aligned,
        sequence_type,
        errors,
        error
    ) = prepare_alignment_and_type(
        files
    )

    if error:

        return jsonify({

            "success": False,

            "message": error

        }), 400

    # Build the hybrid topology
    nj_newick, nj_error = (
        build_neighbor_joining_tree(
            aligned
        )
    )

    if nj_error:

        return jsonify({

            "success": False,

            "message": nj_error

        }), 400

    temp_alignment = None

    try:

        handle = tempfile.NamedTemporaryFile(

            mode="w",

            suffix=".fasta",

            delete=False,

            encoding="utf-8"
        )

        temp_alignment = handle.name

        handle.close()

        write_alignment_file(
            aligned,
            temp_alignment
        )

        hybrid_newick, ml_error = (
            build_ml_tree(

                temp_alignment,

                sequence_type,

                starting_tree=nj_newick
            )
        )

        if ml_error:

            return jsonify({

                "success": False,

                "message": ml_error

            }), 500

        original_clades = (
            get_tree_clades(
                hybrid_newick
            )
        )

        if not original_clades:

            return jsonify({

                "success": False,

                "message":
                    "No internal branches found."

            }), 400

        replicate_clades = []

        support = run_bootstrap_batch(

            aligned,

            original_clades,

            replicate_clades,

            initial
        )

        summary = confidence_summary(
            support
        )

        branches = (
            branch_results_from_support(
                support
            )
        )

        # Adaptive decision
        needs_round2 = (

            any(
                b["support"] < 80
                for b in branches
            )

            and

            initial < maximum
        )

        token = uuid.uuid4().hex

        BOOTSTRAP_JOBS[token] = {

            "aligned":
                aligned,

            "original_newick":
                hybrid_newick,

            "original_clades":
                original_clades,

            "replicate_clades":
                replicate_clades,

            "total":
                initial,

            "batch":
                batch,

            "maximum":
                maximum,

            "stability_threshold":
                threshold,

            "previous_support":
                support
        }

        return jsonify({

            "success":
                True,

            "round":
                1,

            "token":
                token,

            "replicates":
                initial,

            "average_confidence":
                summary[
                    "average_confidence"
                ],

            "branches":
                branches,

            "stable":
                not needs_round2,

            "needs_round2":
                needs_round2,

            "unreliable_regions":
                [
                    b
                    for b in branches
                    if b["support"] < 50
                ],

            "high_confidence":
                summary[
                    "high_confidence"
                ],

            "moderate_confidence":
                summary[
                    "moderate_confidence"
                ],

            "low_confidence":
                summary[
                    "low_confidence"
                ],

            "hybrid_newick":
                hybrid_newick,

            "execution_time":
                round(
                    time.time() - start,
                    4
                ),

            "number_of_sequences":
                len(aligned),

            "alignment_length":
                len(
                    aligned[0]["sequence"]
                )

        })

    finally:

        if (
            temp_alignment
            and
            os.path.exists(
                temp_alignment
            )
        ):

            try:
                os.remove(
                    temp_alignment
                )
            except OSError:
                pass


# ============================================================
# BOOTSTRAP ROUND 2
# ============================================================

@app.route(
    "/bootstrap-round2",
    methods=["POST"]
)
def bootstrap_round2():

    start = time.time()

    token = request.form.get(
        "token",
        ""
    )

    job = BOOTSTRAP_JOBS.get(
        token
    )

    if not job:

        return jsonify({

            "success": False,

            "message":
                "Bootstrap session expired."

        }), 400

    remaining = (
        job["maximum"]
        -
        job["total"]
    )

    if remaining <= 0:

        support = (
            job["previous_support"]
        )

        summary = (
            confidence_summary(
                support
            )
        )

        branches = (
            branch_results_from_support(
                support
            )
        )

        return jsonify({

            "success": True,

            "round": 2,

            "replicates":
                job["total"],

            "branches":
                branches,

            "average_confidence":
                summary[
                    "average_confidence"
                ],

            "stable": False,

            "stopping_reason":
                "Maximum replicate limit reached",

            "unreliable_regions":
                [
                    b
                    for b in branches
                    if b["support"] < 50
                ],

            "hybrid_newick":
                job["original_newick"]

        })

    number = min(
        job["batch"],
        remaining
    )

    support = run_bootstrap_batch(

        job["aligned"],

        job["original_clades"],

        job["replicate_clades"],

        number
    )

    job["total"] += number

    previous = (
        job["previous_support"]
    )

    changes = [

        abs(
            support[c]
            -
            previous.get(
                c,
                support[c]
            )
        )

        for c in support

        if c in previous
    ]

    stable = (

        bool(changes)

        and

        all(
            change
            <= job[
                "stability_threshold"
            ]

            for change in changes
        )
    )

    job["previous_support"] = support

    summary = confidence_summary(
        support
    )

    branches = (
        branch_results_from_support(
            support
        )
    )

    if stable:

        reason = (
            "Branch confidence stabilized"
        )

    elif (
        job["total"]
        >=
        job["maximum"]
    ):

        reason = (
            "Maximum replicate limit reached"
        )

    else:

        reason = (
            "Additional batch completed"
        )

    return jsonify({

        "success":
            True,

        "round":
            2,

        "replicates":
            job["total"],

        "branches":
            branches,

        "average_confidence":
            summary[
                "average_confidence"
            ],

        "stable":
            stable,

        "stopping_reason":
            reason,

        "unreliable_regions":
            [
                b
                for b in branches
                if b["support"] < 50
            ],

        "high_confidence":
            summary[
                "high_confidence"
            ],

        "moderate_confidence":
            summary[
                "moderate_confidence"
            ],

        "low_confidence":
            summary[
                "low_confidence"
            ],

        "hybrid_newick":
            job["original_newick"],

        "execution_time":
            round(
                time.time() - start,
                4
            )
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