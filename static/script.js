// ============================================================
// PHYLOBUILDER
// MULTIPLE FASTA + ALIGNMENT + COLOUR DISPLAY
// ============================================================


let validated = false;


// ============================================================
// SCROLL
// ============================================================

function scrollToUpload() {

    const section =
        document.getElementById(
            "upload-section"
        );

    if (section) {

        section.scrollIntoView({
            behavior: "smooth"
        });

    }
}


// ============================================================
// GET FILES
// ============================================================

function getFiles() {

    const input =
        document.getElementById(
            "fastaFile"
        );

    if (!input) {

        return [];

    }

    return Array.from(
        input.files
    );
}


// ============================================================
// SHOW SELECTED FILES
// ============================================================

function showFileName() {

    const files = getFiles();

    const fileName =
        document.getElementById(
            "fileName"
        );

    const selectedFiles =
        document.getElementById(
            "selectedFiles"
        );

    const alignBtn =
        document.getElementById(
            "alignBtn"
        );

    const njBtn =
        document.getElementById(
            "njBtn"
        );


    validated = false;


    if (alignBtn) {

        alignBtn.disabled = true;

    }


    if (njBtn) {

        njBtn.disabled = true;

    }


    if (files.length === 0) {

        fileName.textContent =
            "No files selected";

        selectedFiles.innerHTML = "";

        return;

    }


    fileName.textContent =
        files.length +
        " FASTA file" +
        (files.length > 1 ? "s" : "") +
        " selected";


    selectedFiles.innerHTML =
        files.map(
            function(file, index) {

                return `
                    <div class="file-chip">
                        ${index + 1}.
                        ${escapeHTML(file.name)}
                    </div>
                `;

            }
        ).join("");

}


// ============================================================
// FORM DATA
// ============================================================

function createFormData() {

    const formData =
        new FormData();

    const files =
        getFiles();


    files.forEach(
        function(file) {

            formData.append(
                "files",
                file
            );

        }
    );


    return formData;
}


// ============================================================
// VALIDATE
// ============================================================

async function validateFiles() {

    const files =
        getFiles();

    const status =
        document.getElementById(
            "validationStatus"
        );


    if (files.length === 0) {

        status.innerHTML = `
            <div class="error-box">
                <h3>Validation Failed</h3>
                <p>
                    Please select at least one FASTA file.
                </p>
            </div>
        `;

        return;

    }


    status.innerHTML = `
        <div class="loading">
            Validating FASTA files...
        </div>
    `;


    try {

        const response =
            await fetch(
                "/validate",
                {
                    method: "POST",
                    body: createFormData()
                }
            );


        const data =
            await response.json();


        if (
            !response.ok ||
            !data.success
        ) {

            const message =
                data.message ||
                (
                    data.errors
                        ? data.errors.join("<br>")
                        : "Validation failed."
                );


            throw new Error(
                message
            );

        }


        const stats =
            data.statistics;


        // IMPORTANT:
        // The Flask API returns
        // statistics.number_of_sequences
        // NOT data.count.


        document.getElementById(
            "sequenceCount"
        ).textContent =
            stats.number_of_sequences;


        document.getElementById(
            "sequenceType"
        ).textContent =
            stats.sequence_type;


        document.getElementById(
            "totalLength"
        ).textContent =
            stats.total_length;


        document.getElementById(
            "fileCount"
        ).textContent =
            stats.number_of_files;


        document.getElementById(
            "summarySequences"
        ).textContent =
            stats.number_of_sequences;


        status.innerHTML = `
            <div class="success-box">

                <h3>
                    ✓ FASTA Validation Successful
                </h3>

                <p>
                    ${stats.number_of_sequences}
                    valid sequences detected.
                </p>

                <p>
                    Files:
                    ${stats.number_of_files}
                </p>

                ${
                    data.warnings &&
                    data.warnings.length
                    ?
                    `
                    <small>
                        ${escapeHTML(
                            data.warnings.join(" | ")
                        )}
                    </small>
                    `
                    :
                    ""
                }

            </div>
        `;


        validated = true;


        document.getElementById(
            "alignBtn"
        ).disabled = false;


    }
    catch (error) {

        console.error(
            "Validation error:",
            error
        );


        status.innerHTML = `
            <div class="error-box">

                <h3>
                    Validation Failed
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


// ============================================================
// RUN ALIGNMENT
// ============================================================

async function runAlignment() {

    const files =
        getFiles();


    const preview =
        document.getElementById(
            "alignmentPreview"
        );


    const resultBox =
        document.getElementById(
            "alignmentResult"
        );


    const meta =
        document.getElementById(
            "alignmentMeta"
        );


    const button =
        document.getElementById(
            "alignBtn"
        );


    if (files.length === 0) {

        resultBox.innerHTML = `
            <div class="error-box">

                <h3>
                    Alignment Failed
                </h3>

                <p>
                    Please select FASTA files first.
                </p>

            </div>
        `;

        return;

    }


    button.disabled = true;


    preview.innerHTML = `
        <div class="loading">
            Running Multiple Sequence Alignment...
        </div>
    `;


    try {

        const response =
            await fetch(
                "/align",
                {
                    method: "POST",
                    body: createFormData()
                }
            );


        const data =
            await response.json();


        if (
            !response.ok ||
            !data.success
        ) {

            const message =
                data.message ||
                (
                    data.errors
                        ? data.errors.join("<br>")
                        : "Alignment failed."
                );


            throw new Error(
                message
            );

        }


        // ====================================================
        // ALIGNMENT IS AN ARRAY
        //
        // Example:
        //
        // [
        //   {
        //      name: "Human",
        //      sequence: "ATGC..."
        //   }
        // ]
        // ====================================================


        const alignment =
            data.alignment;


        if (
            !alignment ||
            !Array.isArray(alignment)
        ) {

            throw new Error(
                "No alignment data returned by the server."
            );

        }


        renderAlignment(
            alignment
        );


        meta.textContent =
            data.number_of_sequences +
            " sequences | " +
            data.alignment_length +
            " positions | " +
            data.execution_time +
            " seconds";


        resultBox.innerHTML = `
            <div class="result-card">

                <h3>
                    ✓ Multiple Sequence Alignment Completed
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


        document.getElementById(
            "summarySequences"
        ).textContent =
            data.number_of_sequences;


        document.getElementById(
            "summaryAlignment"
        ).textContent =
            data.alignment_length +
            " positions";


        const njButton =
            document.getElementById(
                "njBtn"
            );


        if (data.number_of_sequences >= 3) {

            njButton.disabled = false;

        }


        document.getElementById(
            "alignment-section"
        ).scrollIntoView({
            behavior: "smooth"
        });


    }
    catch (error) {

        console.error(
            "Alignment error:",
            error
        );


        preview.innerHTML = `
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


        resultBox.innerHTML = "";

    }


    button.disabled = !validated;

}


// ============================================================
// COLOUR-CODED ALIGNMENT
// ============================================================

function renderAlignment(
    rows
) {

    const preview =
        document.getElementById(
            "alignmentPreview"
        );


    if (
        !rows ||
        rows.length === 0
    ) {

        preview.innerHTML =
            "<p>No alignment available.</p>";

        return;

    }


    const alignmentLength =
        rows[0].sequence.length;


    // ========================================================
    // CREATE CONSENSUS
    // ========================================================

    let consensus = "";


    for (
        let position = 0;
        position < alignmentLength;
        position++
    ) {

        const counts = {};


        rows.forEach(
            function(row) {

                const character =
                    (
                        row.sequence[position] ||
                        "-"
                    ).toUpperCase();


                if (
                    character !== "-"
                ) {

                    counts[character] =
                        (
                            counts[character] ||
                            0
                        ) + 1;

                }

            }
        );


        const sorted =
            Object.entries(
                counts
            ).sort(
                function(a, b) {

                    return b[1] - a[1];

                }
            );


        if (sorted.length > 0) {

            consensus +=
                sorted[0][0];

        }
        else {

            consensus += "-";

        }

    }


    // ========================================================
    // CONSENSUS ROW
    // ========================================================

    let consensusHTML = "";


    for (
        let i = 0;
        i < consensus.length;
        i++
    ) {

        const character =
            consensus[i];


        consensusHTML += `
            <span class="consensus-char">
                ${escapeHTML(character)}
            </span>
        `;

    }


    // ========================================================
    // SEQUENCE ROWS
    // ========================================================

    let rowsHTML = "";


    rows.forEach(
        function(row) {

            let sequenceHTML = "";


            for (
                let i = 0;
                i < row.sequence.length;
                i++
            ) {

                const character =
                    row.sequence[i]
                    .toUpperCase();


                const isConserved =
                    character !== "-" &&
                    character === consensus[i];


                const className =
                    getBaseClass(
                        character,
                        isConserved
                    );


                sequenceHTML += `
                    <span
                        class="base ${className}"
                        title="Position ${i + 1}"
                    >
                        ${escapeHTML(character)}
                    </span>
                `;

            }


            rowsHTML += `
                <div class="alignment-row">

                    <div
                        class="sequence-name"
                        title="${escapeHTML(row.name)}"
                    >
                        ${escapeHTML(row.name)}
                    </div>

                    <div class="sequence-track">
                        ${sequenceHTML}
                    </div>

                </div>
            `;

        }
    );


    // ========================================================
    // DISPLAY
    // ========================================================

    preview.innerHTML = `

        <div class="alignment-container">

            <div class="alignment-row consensus-row">

                <div class="sequence-name">
                    CONSENSUS
                </div>

                <div class="sequence-track">
                    ${consensusHTML}
                </div>

            </div>

            ${rowsHTML}

        </div>

    `;

}


// ============================================================
// BASE COLOUR
// ============================================================

function getBaseClass(
    character,
    conserved
) {

    let className = "";


    if (character === "A") {

        className = "base-A";

    }
    else if (character === "T") {

        className = "base-T";

    }
    else if (character === "G") {

        className = "base-G";

    }
    else if (character === "C") {

        className = "base-C";

    }
    else if (character === "U") {

        className = "base-U";

    }
    else if (character === "-") {

        className = "base-gap";

    }
    else {

        className = "base-other";

    }


    if (conserved) {

        className += " conserved";

    }


    return className;

}


// ============================================================
// NEIGHBOR JOINING TREE
// ============================================================

async function buildNJTree() {

    const files =
        getFiles();


    const resultBox =
        document.getElementById(
            "njResult"
        );


    const button =
        document.getElementById(
            "njBtn"
        );


    if (files.length === 0) {

        resultBox.innerHTML = `
            <div class="error-box">

                <h3>
                    Tree Construction Failed
                </h3>

                <p>
                    Please upload FASTA files first.
                </p>

            </div>
        `;

        return;

    }


    button.disabled = true;


    resultBox.innerHTML = `
        <div class="loading">
            Building Neighbor Joining Tree...
        </div>
    `;


    try {

        const response =
            await fetch(
                "/build-tree",
                {
                    method: "POST",
                    body: createFormData()
                }
            );


        const data =
            await response.json();


        if (
            !response.ok ||
            !data.success
        ) {

            const message =
                data.message ||
                (
                    data.errors
                        ? data.errors.join("<br>")
                        : "Tree construction failed."
                );


            throw new Error(
                message
            );

        }


        resultBox.innerHTML = `

            <div class="result-card">

                <h3>
                    ✓ Neighbor Joining Tree Generated
                </h3>

                <p>
                    <strong>
                        Method:
                    </strong>

                    ${data.method}
                </p>

                <p>
                    <strong>
                        Sequences:
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

                <pre class="newick-box">
${escapeHTML(data.newick)}
                </pre>

            </div>

        `;


        document.getElementById(
            "summaryTree"
        ).textContent =
            "NJ Ready";


    }
    catch (error) {

        console.error(
            "Tree error:",
            error
        );


        resultBox.innerHTML = `
            <div class="error-box">

                <h3>
                    Tree Construction Error
                </h3>

                <p>
                    ${escapeHTML(
                        error.message
                    )}
                </p>

            </div>
        `;

    }


    button.disabled = false;

}


// ============================================================
// ESCAPE HTML
// ============================================================

function escapeHTML(
    value
) {

    const div =
        document.createElement(
            "div"
        );


    div.textContent =
        String(
            value ?? ""
        );


    return div.innerHTML;

}