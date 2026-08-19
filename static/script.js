function scrollToUpload() {
    document.getElementById("upload-section")
        .scrollIntoView({ behavior: "smooth" });
}


function showFileName() {

    const fileInput = document.getElementById("fastaFile");
    const fileName = document.getElementById("fileName");

    if (fileInput.files.length > 0) {
        fileName.textContent =
            "Selected: " + fileInput.files[0].name;
    }
}


async function validateFile() {

    const fileInput = document.getElementById("fastaFile");

    if (fileInput.files.length === 0) {
        alert("Please select a FASTA file.");
        return;
    }

    const formData = new FormData();

    formData.append("file", fileInput.files[0]);

    try {

        const response = await fetch("/validate", {
            method: "POST",
            body: formData
        });

        const data = await response.json();

        if (!data.success) {

            alert(
                data.errors
                    ? data.errors.join("\n")
                    : data.message
            );

            return;
        }

        displayResults(data);

    } catch (error) {

        alert("Unable to connect to the server.");

        console.error(error);
    }
}


function displayResults(data) {

    document.getElementById("results")
        .classList.remove("hidden");

    document.getElementById("sequenceCount")
        .textContent = data.count;

    const table =
        document.getElementById("sequenceTable");

    table.innerHTML = "";

    let totalLength = 0;

    data.sequences.forEach((seq, index) => {

        totalLength += seq.length;

        const row = document.createElement("tr");

        row.innerHTML = `
            <td>${index + 1}</td>
            <td>${seq.name}</td>
            <td>${seq.length}</td>
        `;

        table.appendChild(row);
    });

    document.getElementById("totalLength")
        .textContent = totalLength;

    detectSequenceType(data.sequences);
}


function detectSequenceType(sequences) {

    const dnaChars = /^[ACGTN]+$/;
    const rnaChars = /^[ACGUN]+$/;

    let type = "Protein";

    const allDNA = sequences.every(
        seq => dnaChars.test(seq.sequence)
    );

    const allRNA = sequences.every(
        seq => rnaChars.test(seq.sequence)
    );

    if (allDNA) {
        type = "DNA";
    } else if (allRNA) {
        type = "RNA";
    }

    document.getElementById("sequenceType")
        .textContent = type;
}