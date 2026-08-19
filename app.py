from flask import Flask, render_template, request, jsonify
from Bio import SeqIO
from io import StringIO

app = Flask(__name__)


# ---------------------------------------------------------
# FASTA VALIDATION
# ---------------------------------------------------------

def validate_fasta(text):
    sequences = []
    errors = []

    try:
        records = list(SeqIO.parse(StringIO(text), "fasta"))

        # No sequences found
        if len(records) == 0:
            errors.append("No FASTA sequences found.")
            return sequences, errors

        # Check every sequence
        for record in records:

            sequence = str(record.seq).upper()

            # Check for empty sequence
            if len(sequence) == 0:
                error_message = str(record.id) + ": Empty sequence."
                errors.append(error_message)
                continue

            # Store sequence information
            sequence_data = {
                "name": record.id,
                "length": len(sequence),
                "sequence": sequence
            }

            sequences.append(sequence_data)

    except Exception as error:

        error_message = (
            "Invalid FASTA format: " + str(error)
        )

        errors.append(error_message)

    return sequences, errors


# ---------------------------------------------------------
# HOME PAGE
# ---------------------------------------------------------

@app.route("/")
def home():
    return render_template("index.html")


# ---------------------------------------------------------
# VALIDATE FASTA FILE
# ---------------------------------------------------------

@app.route("/validate", methods=["POST"])
def validate():

    # Check whether a file was uploaded
    if "file" not in request.files:

        return jsonify({
            "success": False,
            "message": "Please upload a FASTA file."
        })

    file = request.files["file"]

    # Check whether the user selected a file
    if file.filename == "":

        return jsonify({
            "success": False,
            "message": "No file selected."
        })

    try:

        # Read uploaded file
        file_content = file.read()

        # Convert bytes to text
        text = file_content.decode("utf-8")

    except UnicodeDecodeError:

        return jsonify({
            "success": False,
            "message": "The uploaded file is not a valid UTF-8 text file."
        })

    # Validate FASTA
    sequences, errors = validate_fasta(text)

    # If validation failed
    if len(errors) > 0:

        return jsonify({
            "success": False,
            "errors": errors
        })

    # Calculate total length
    total_length = 0

    for sequence in sequences:
        total_length += sequence["length"]

    # Return successful result
    return jsonify({
        "success": True,
        "count": len(sequences),
        "total_length": total_length,
        "sequences": sequences,
        "message": "FASTA validation successful."
    })


# ---------------------------------------------------------
# RUN FLASK APPLICATION
# ---------------------------------------------------------

if __name__ == "__main__":
    app.run(
        host="127.0.0.1",
        port=5000,
        debug=True
    )