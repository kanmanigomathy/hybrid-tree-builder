// ============================================================
// PHYLOBUIDER - MAIN JAVASCRIPT
// FASTA Validation → Alignment → Neighbor Joining
// ============================================================


// ============================================================
// SCROLL TO UPLOAD
// ============================================================

function scrollToUpload() {

    const section = document.getElementById("upload-section");

    if (section) {

        section.scrollIntoView({
            behavior: "smooth"
        });

    }
}


// ============================================================
// SHOW SELECTED FILE NAME
// ============================================================

function showFileName() {

    const fileInput =
        document.getElementById("fastaFile");

    const fileName =
        document.getElementById("fileName");

    if (!fileInput || !fileName) {
        return;
    }

    if (fileInput.files.length > 0) {

        fileName.textContent =
            "Selected: " +
            fileInput.files[0].name;

    } else {

        fileName.textContent =
            "No file selected";

    }
}


// ============================================================
// FASTA VALIDATION
// ============================================================

async function validateFile() {

    const fileInput =
        document.getElementById("fastaFile");

    if (!fileInput || fileInput.files.length === 0) {

        alert("Please select a FASTA file.");

        return;
    }


    const formData =
        new FormData();

    formData.append(
        "file",
        fileInput.files[0]
    );


    try {

        const response =
            await fetch("/validate", {

                method: "POST",

                body: formData

            });


        const data =
            await response.json();


        if (!data.success) {

            alert(
                data.errors
                    ? data.errors.join("\n")
                    : data.message || "Validation failed."
            );

            return;
        }


        displayResults(data);


        // Enable alignment button
        const alignBtn =
            document.getElementById("alignBtn");

        if (alignBtn) {
            alignBtn.disabled = false;
        }

    }

    catch (error) {

        console.error(
            "Validation error:",
            error
        );

        alert(
            "Unable to connect to the Flask server."
        );

    }
}


// ============================================================
// DISPLAY VALIDATION RESULTS
// ============================================================

function displayResults(data) {

    const results =
        document.getElementById("results");

    const sequenceCount =
        document.getElementById("sequenceCount");

    const table =
        document.getElementById("sequenceTable");

    const totalLength =
        document.getElementById("totalLength");


    if (results) {

        results.classList.remove("hidden");

    }


    if (sequenceCount) {

        sequenceCount.textContent =
            data.count;

    }


    if (!table) {
        return;
    }


    table.innerHTML = "";


    let total =
        0;


    data.sequences.forEach(
        (seq, index) => {

            total +=
                seq.length;


            const row =
                document.createElement("tr");


            row.innerHTML = `

                <td>
                    ${index + 1}
                </td>

                <td>
                    ${escapeHTML(seq.name)}
                </td>

                <td>
                    ${seq.length}
                </td>

            `;


            table.appendChild(row);

        }
    );


    if (totalLength) {

        totalLength.textContent =
            total;

    }


    detectSequenceType(
        data.sequences
    );
}


// ============================================================
// DETECT DNA / RNA / PROTEIN
// ============================================================

function detectSequenceType(sequences) {

    const dnaChars =
        /^[ACGTN]+$/i;

    const rnaChars =
        /^[ACGUN]+$/i;


    let type =
        "Protein";


    const allDNA =
        sequences.every(
            seq =>
                dnaChars.test(
                    seq.sequence
                )
        );


    const allRNA =
        sequences.every(
            seq =>
                rnaChars.test(
                    seq.sequence
                )
        );


    if (allDNA) {

        type =
            "DNA";

    }

    else if (allRNA) {

        type =
            "RNA";

    }


    const typeBox =
        document.getElementById(
            "sequenceType"
        );


    if (typeBox) {

        typeBox.textContent =
            type;

    }
}


// ============================================================
// MULTIPLE SEQUENCE ALIGNMENT
// ============================================================

async function runAlignment() {

    const fileInput =
        document.getElementById(
            "fastaFile"
        );


    const resultBox =
        document.getElementById(
            "alignmentResult"
        );


    const previewBox =
        document.getElementById(
            "alignmentPreview"
        );


    const alignBtn =
        document.getElementById(
            "alignBtn"
        );


    // --------------------------------------------------------
    // CHECK FILE
    // --------------------------------------------------------

    if (
        !fileInput ||
        fileInput.files.length === 0
    ) {

        alert(
            "Please select a FASTA file first."
        );

        return;
    }


    // --------------------------------------------------------
    // DISABLE BUTTON
    // --------------------------------------------------------

    if (alignBtn) {

        alignBtn.disabled =
            true;

    }


    // --------------------------------------------------------
    // SHOW LOADING
    // --------------------------------------------------------

    if (resultBox) {

        resultBox.innerHTML = `

            <div class="loading">

                <h3>
                    Running Multiple Sequence Alignment...
                </h3>

                <p>
                    Please wait.
                </p>

            </div>

        `;

    }


    if (previewBox) {

        previewBox.innerHTML = `

            <div class="loading">

                Running alignment...

            </div>

        `;

    }


    // --------------------------------------------------------
    // CREATE FORM DATA
    // --------------------------------------------------------

    const formData =
        new FormData();


    formData.append(
        "file",
        fileInput.files[0]
    );


    try {

        console.log(
            "Sending FASTA file:",
            fileInput.files[0].name
        );


        // ----------------------------------------------------
        // SEND FILE TO FLASK
        // ----------------------------------------------------

        const response =
            await fetch(
                "/align",
                {

                    method: "POST",

                    body: formData

                }
            );


        console.log(
            "Alignment HTTP status:",
            response.status
        );


        // ----------------------------------------------------
        // READ RESPONSE
        // ----------------------------------------------------

        const responseText =
            await response.text();


        console.log(
            "Alignment server response:",
            responseText
        );


        let data;


        try {

            data =
                JSON.parse(
                    responseText
                );

        }

        catch (jsonError) {

            throw new Error(
                "Flask returned an invalid response: " +
                responseText
            );

        }


        // ----------------------------------------------------
        // CHECK ALIGNMENT SUCCESS
        // ----------------------------------------------------

        if (!data.success) {

            const message =
                data.message ||
                data.error ||
                "Alignment failed.";


            if (resultBox) {

                resultBox.innerHTML = `

                    <div class="error-box">

                        <h3>
                            Alignment Failed
                        </h3>

                        <p>
                            ${escapeHTML(message)}
                        </p>

                    </div>

                `;

            }


            if (previewBox) {

                previewBox.innerHTML = `

                    <div class="error-box">

                        <strong>
                            Alignment Failed
                        </strong>

                        <p>
                            ${escapeHTML(message)}
                        </p>

                    </div>

                `;

            }


            return;
        }


        // ----------------------------------------------------
        // BUILD ALIGNMENT DISPLAY
        // ----------------------------------------------------

        let alignmentHTML = "";

if (data.alignment) {

    // Flask returns alignment as an array
    if (Array.isArray(data.alignment)) {

        data.alignment.forEach(seq => {

            alignmentHTML += `
                <div class="alignment-row">
                    <strong>${escapeHTML(seq.name)}</strong>
                    &nbsp;&nbsp;
                    ${escapeHTML(seq.sequence)}
                </div>
            `;

        });

    }

    // Backup in case Flask returns a string
    else if (typeof data.alignment === "string") {

        const lines = data.alignment
            .split("\n")
            .filter(line => line.trim() !== "");

        lines.forEach(line => {

            alignmentHTML += `
                <div class="alignment-row">
                    ${escapeHTML(line)}
                </div>
            `;

        });

    }

}


        // ----------------------------------------------------
        // SHOW ALIGNMENT PREVIEW
        // ----------------------------------------------------

        if (previewBox) {

            previewBox.innerHTML = `

                <div class="success-box">

                    <h3>
                        ✓ Multiple Sequence Alignment Completed
                    </h3>


                    <div class="alignment-summary">

                        <p>
                            <strong>
                                Number of sequences:
                            </strong>

                            ${data.number_of_sequences}
                        </p>


                        <p>
                            <strong>
                                Alignment length:
                            </strong>

                            ${data.alignment_length}
                        </p>


                        <p>
                            <strong>
                                Execution time:
                            </strong>

                            ${data.execution_time}
                            seconds
                        </p>

                    </div>


                    <h4>
                        Alignment
                    </h4>


                    <div class="alignment-box">

                        ${alignmentHTML}

                    </div>

                </div>

            `;

        }


        // ----------------------------------------------------
        // SHOW RESULT CARD
        // ----------------------------------------------------

        if (resultBox) {

            resultBox.innerHTML = `

                <div class="result-card">

                    <h3>
                        ✓ Alignment Completed Successfully
                    </h3>


                    <p>

                        <strong>
                            Number of sequences:
                        </strong>

                        ${data.number_of_sequences}

                    </p>


                    <p>

                        <strong>
                            Alignment length:
                        </strong>

                        ${data.alignment_length}

                    </p>


                    <p>

                        <strong>
                            Execution time:
                        </strong>

                        ${data.execution_time}
                        seconds

                    </p>

                </div>

            `;

        }


        // ----------------------------------------------------
        // ENABLE NEIGHBOR JOINING BUTTON
        // ----------------------------------------------------

        const njBtn =
            document.getElementById(
                "njBtn"
            );


        if (njBtn) {

            njBtn.disabled =
                false;

        }

    }


    catch (error) {

        console.error(
            "REAL ALIGNMENT ERROR:",
            error
        );


        // ----------------------------------------------------
        // SHOW ACTUAL ERROR
        // ----------------------------------------------------

        if (resultBox) {

            resultBox.innerHTML = `

                <div class="error-box">

                    <h3>
                        Alignment Error
                    </h3>

                    <p>
                        ${escapeHTML(
                            error.message
                        )}
                    </p>

                </div>

            `;

        }


        if (previewBox) {

            previewBox.innerHTML = `

                <div class="error-box">

                    <h3>
                        Alignment Error
                    </h3>

                    <p>
                        ${escapeHTML(
                            error.message
                        )}
                    </p>

                </div>

            `;

        }

    }


    finally {

        if (alignBtn) {

            alignBtn.disabled =
                false;

        }

    }
}


// ============================================================
// ESCAPE HTML
// ============================================================

function escapeHTML(text) {

    const div =
        document.createElement(
            "div"
        );


    div.textContent =
        text;


    return div.innerHTML;
}


// ============================================================
// NEIGHBOR JOINING TREE
// ============================================================

async function buildNJTree() {

    const fileInput =
        document.getElementById(
            "fastaFile"
        );


    const resultBox =
        document.getElementById(
            "njResult"
        );


    const njBtn =
        document.getElementById(
            "njBtn"
        );


    // --------------------------------------------------------
    // CHECK FILE
    // --------------------------------------------------------

    if (
        !fileInput ||
        fileInput.files.length === 0
    ) {

        alert(
            "Please select a FASTA file first."
        );

        return;
    }


    // --------------------------------------------------------
    // DISABLE NJ BUTTON
    // --------------------------------------------------------

    if (njBtn) {

        njBtn.disabled =
            true;

    }


    // --------------------------------------------------------
    // SHOW LOADING
    // --------------------------------------------------------

    if (resultBox) {

        resultBox.innerHTML = `

            <div class="loading">

                <h3>
                    Building Neighbor Joining Tree...
                </h3>

                <p>
                    Please wait.
                </p>

            </div>

        `;

    }


    // --------------------------------------------------------
    // CREATE FORM DATA
    // --------------------------------------------------------

    const formData =
        new FormData();


    formData.append(
        "file",
        fileInput.files[0]
    );


    try {

        const response =
            await fetch(
                "/build-tree",
                {

                    method: "POST",

                    body: formData

                }
            );


        const responseText =
            await response.text();


        console.log(
            "NJ server response:",
            responseText
        );


        let data;


        try {

            data =
                JSON.parse(
                    responseText
                );

        }

        catch (error) {

            throw new Error(
                "Invalid response from Flask."
            );

        }


        // ----------------------------------------------------
        // TREE FAILED
        // ----------------------------------------------------

        if (!data.success) {

            const message =
                data.message ||
                data.error ||
                "Tree construction failed.";


            if (resultBox) {

                resultBox.innerHTML = `

                    <div class="error-box">

                        <h3>
                            Tree Construction Failed
                        </h3>

                        <p>
                            ${escapeHTML(message)}
                        </p>

                    </div>

                `;

            }


            return;
        }


        // ----------------------------------------------------
        // SHOW NJ RESULT
        // ----------------------------------------------------

        if (resultBox) {

            resultBox.innerHTML = `

                <div class="result-card">

                    <h3>
                        ✓ Neighbor Joining Tree Generated
                    </h3>


                    <p>

                        <strong>
                            Method:
                        </strong>

                        ${escapeHTML(
                            data.method ||
                            "Neighbor Joining"
                        )}

                    </p>


                    <p>

                        <strong>
                            Number of sequences:
                        </strong>

                        ${data.number_of_sequences}

                    </p>


                    <p>

                        <strong>
                            Execution time:
                        </strong>

                        ${data.execution_time}
                        seconds

                    </p>


                    <h4>
                        Newick Representation
                    </h4>


                    <pre class="newick-box">${escapeHTML(
                        data.newick || ""
                    )}</pre>

                </div>

            `;

        }

    }


    catch (error) {

        console.error(
            "NJ error:",
            error
        );


        if (resultBox) {

            resultBox.innerHTML = `

                <div class="error-box">

                    <h3>
                        Neighbor Joining Error
                    </h3>

                    <p>
                        ${escapeHTML(
                            error.message
                        )}
                    </p>

                </div>

            `;

        }

    }


    finally {

        if (njBtn) {

            njBtn.disabled =
                false;

        }

    }
}