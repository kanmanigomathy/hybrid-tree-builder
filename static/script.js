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

// ============================================================
// PHASE 2 — HYBRID TREE
// NEIGHBOR JOINING + MAXIMUM LIKELIHOOD
// ============================================================

async function buildNJTree() {

    const files = getFiles();

    const resultBox =
        document.getElementById("njResult");

    const button =
        document.getElementById("njBtn");


    // --------------------------------------------------------
    // CHECK FILES
    // --------------------------------------------------------

    if (!files || files.length === 0) {

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


    // --------------------------------------------------------
    // DISABLE BUTTON
    // --------------------------------------------------------

    button.disabled = true;


    resultBox.innerHTML = `
        <div class="loading">

            <p>
                Building Hybrid Phylogenetic Tree...
            </p>

            <p>
                Step 1: Neighbor Joining
            </p>

            <p>
                Step 2: Maximum Likelihood refinement
            </p>

        </div>
    `;


    try {

        // ----------------------------------------------------
        // SEND FASTA FILES TO FLASK
        // ----------------------------------------------------

        const response =
            await fetch(
                "/build-tree",
                {
                    method: "POST",
                    body: createFormData()
                }
            );


        // ----------------------------------------------------
        // CHECK SERVER RESPONSE
        // ----------------------------------------------------

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
                        : "Hybrid tree construction failed."
                );


            throw new Error(message);
        }


        // ====================================================
        // GET TREE DATA
        // ====================================================

        /*
         * The backend may use different names depending on
         * the current app.py version.
         *
         * Therefore we support several possible keys.
         */

        const njNewick =
            data.nj_newick ||
            data.nj_tree ||
            data.newick ||
            "";


        const mlNewick =
            data.ml_newick ||
            data.ml_tree ||
            data.maximum_likelihood ||
            "";


        const finalNewick =
            data.final_newick ||
            data.final_tree ||
            data.hybrid_newick ||
            data.hybrid_tree ||
            mlNewick ||
            njNewick ||
            "";


        // ----------------------------------------------------
        // EXECUTION TIME
        // ----------------------------------------------------

        const executionTime =
            data.execution_time !== undefined
                ? data.execution_time
                : "N/A";


        // ----------------------------------------------------
        // NUMBER OF SEQUENCES
        // ----------------------------------------------------

        const numberOfSequences =
            data.number_of_sequences ||
            files.length;


        // ====================================================
        // BUILD NJ DISPLAY
        // ====================================================

        let njHTML = "";

        if (njNewick) {

            njHTML = `
                <div class="tree-stage">

                    <h4>
                        🌿 Step 1 — Neighbor Joining Tree
                    </h4>

                    <p class="tree-description">
                        Fast distance-based initial phylogenetic tree.
                    </p>

                    <pre class="newick-box">${escapeHTML(
                        njNewick
                    )}</pre>

                </div>
            `;

        }
        else {

            njHTML = `
                <div class="tree-stage">

                    <h4>
                        🌿 Step 1 — Neighbor Joining Tree
                    </h4>

                    <p>
                        NJ tree data was not returned by the server.
                    </p>

                </div>
            `;

        }


        // ====================================================
        // BUILD ML DISPLAY
        // ====================================================

        let mlHTML = "";

        if (mlNewick) {

            mlHTML = `
                <div class="tree-stage">

                    <h4>
                        🧬 Step 2 — Maximum Likelihood Tree
                    </h4>

                    <p class="tree-description">
                        The Neighbor Joining topology is evaluated
                        and refined using Maximum Likelihood.
                    </p>

                    <pre class="newick-box">${escapeHTML(
                        mlNewick
                    )}</pre>

                </div>
            `;

        }
        else {

            mlHTML = `
                <div class="tree-stage">

                    <h4>
                        🧬 Step 2 — Maximum Likelihood Tree
                    </h4>

                    <div class="warning-box">

                        <strong>
                            ML tree not returned
                        </strong>

                        <p>
                            The FastTree Maximum Likelihood
                            result is not present in the server response.
                        </p>

                    </div>

                </div>
            `;

        }


        // ====================================================
        // BUILD FINAL HYBRID TREE DISPLAY
        // ====================================================

        let finalHTML = "";

        if (finalNewick) {

            finalHTML = `
                <div class="tree-stage final-tree-stage">

                    <h4>
                        🌳 Final Hybrid Phylogenetic Tree
                    </h4>

                    <p class="tree-description">
                        Final tree obtained after combining
                        Neighbor Joining initialization with
                        Maximum Likelihood refinement.
                    </p>

                    <pre class="newick-box final-newick">${escapeHTML(
                        finalNewick
                    )}</pre>

                </div>
            `;

        }


        // ====================================================
        // COMPLETE RESULT
        // ====================================================

        resultBox.innerHTML = `

            <div class="result-card">

                <h3>
                    ✓ Hybrid Phylogenetic Tree Generated
                </h3>


                <p>

                    <strong>
                        Method:
                    </strong>

                    ${
                        data.method ||
                        "Hybrid Neighbor Joining + Maximum Likelihood"
                    }

                </p>


                <p>

                    <strong>
                        Sequences:
                    </strong>

                    ${numberOfSequences}

                </p>


                <p>

                    <strong>
                        Execution time:
                    </strong>

                    ${executionTime}

                    seconds

                </p>


                <!-- ===================================== -->
                <!-- WORKFLOW -->
                <!-- ===================================== -->

                <div class="tree-workflow">

                    <div class="workflow-step active">
                        <span>1</span>
                        <strong>
                            Neighbor Joining
                        </strong>
                        <small>
                            Initial tree
                        </small>
                    </div>


                    <div class="workflow-arrow">
                        ↓
                    </div>


                    <div class="workflow-step active">
                        <span>2</span>
                        <strong>
                            Maximum Likelihood
                        </strong>
                        <small>
                            Tree refinement
                        </small>
                    </div>


                    <div class="workflow-arrow">
                        ↓
                    </div>


                    <div class="workflow-step final">
                        <span>3</span>
                        <strong>
                            Final Hybrid Tree
                        </strong>
                        <small>
                            Final result
                        </small>
                    </div>

                </div>


                <!-- ===================================== -->
                <!-- NJ TREE -->
                <!-- ===================================== -->

                ${njHTML}


                <!-- ===================================== -->
                <!-- ML TREE -->
                <!-- ===================================== -->

                ${mlHTML}


                <!-- ===================================== -->
                <!-- FINAL TREE -->
                <!-- ===================================== -->

                ${finalHTML}

            </div>

        `;


        // ====================================================
        // UPDATE SUMMARY
        // ====================================================

        const summaryTree =
            document.getElementById("summaryTree");


        if (summaryTree) {

            summaryTree.textContent =
                mlNewick
                    ? "Hybrid Tree Ready"
                    : "NJ Tree Ready";

        }


        // ====================================================
        // SCROLL TO TREE SECTION
        // ====================================================

        const treeSection =
            document.getElementById("tree-section");


        if (treeSection) {

            treeSection.scrollIntoView({
                behavior: "smooth"
            });

        }

    }


    catch (error) {

        console.error(
            "Hybrid tree error:",
            error
        );


        resultBox.innerHTML = `

            <div class="error-box">

                <h3>
                    Hybrid Tree Construction Error
                </h3>

                <p>
                    ${escapeHTML(
                        error.message
                    )}
                </p>

            </div>

        `;

    }


    // --------------------------------------------------------
    // ENABLE BUTTON AGAIN
    // --------------------------------------------------------

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