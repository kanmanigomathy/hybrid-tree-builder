let validated = false;
let alignmentReady = false;

let njReady = false;
let mlReady = false;
let hybridReady = false;
let bootstrapReady = false;

let phyloTreeData = null;
let bootstrapData = null;

let bootstrapToken = null;

let treeZoom = 1;


// ============================================================
// BASIC HELPERS
// ============================================================

function $(id) {

    return document.getElementById(id);

}


function scrollToUpload() {

    $("upload-section")?.scrollIntoView({
        behavior: "smooth"
    });

}


function getFiles() {

    const input = $("fastaFile");

    return input
        ? Array.from(input.files)
        : [];

}


function createFormData() {

    const fd = new FormData();

    getFiles().forEach(file => {

        fd.append(
            "files",
            file
        );

    });

    return fd;
}


function escapeHTML(value) {

    const div =
        document.createElement(
            "div"
        );

    div.textContent =
        String(value ?? "");

    return div.innerHTML;

}


function setText(id, value) {

    const el = $(id);

    if (el) {

        el.textContent =
            String(value ?? "");

    }

}


function showError(
    element,
    title,
    message
) {

    if (!element) return;

    element.innerHTML = `

        <div class="error-box">

            <h3>
                ${escapeHTML(title)}
            </h3>

            <p>
                ${escapeHTML(message)}
            </p>

        </div>

    `;

}


function showSuccess(
    element,
    title,
    message
) {

    if (!element) return;

    element.innerHTML = `

        <div class="success-box">

            <h3>
                ${escapeHTML(title)}
            </h3>

            <p>
                ${escapeHTML(message)}
            </p>

        </div>

    `;

}


function getNumberValue(
    id,
    fallback
) {

    const el = $(id);

    const value =
        Number(el?.value);

    return Number.isFinite(value)
        ? value
        : fallback;

}


function clamp(
    value,
    min,
    max
) {

    return Math.max(
        min,
        Math.min(max, value)
    );

}


// ============================================================
// RESET
// ============================================================

function resetWorkflow() {

    validated = false;

    alignmentReady = false;

    njReady = false;

    mlReady = false;

    hybridReady = false;

    bootstrapReady = false;

    phyloTreeData = null;

    bootstrapData = null;

    bootstrapToken = null;


    [
        "alignBtn",
        "njBtn",
        "mlBtn",
        "hybridBtn",
        "bootstrapBtn"
    ].forEach(id => {

        if ($(id)) {

            $(id).disabled = true;

        }

    });


    $("njNewick").textContent =
        "NJ tree not built yet.";

    $("mlNewick").textContent =
        "ML tree not built yet.";

    $("hybridNewick").textContent =
        "Hybrid tree not generated yet.";


    $("treeCanvas").innerHTML = `

        <div class="empty-viz">

            The interactive tree will appear
            after hybrid tree generation.

        </div>

    `;


    $("branchDetails").textContent =
        "Click an internal node to view branch confidence.";


    $("bootstrapRound1").innerHTML = "";

    $("bootstrapRound2").innerHTML = "";

    $("bootstrapResult").innerHTML = "";

    $("finalResults").innerHTML = "";

    setText(
        "summaryTree",
        "—"
    );

    setText(
        "summaryBootstrap",
        "Next"
    );

}


// ============================================================
// FILE SELECTION
// ============================================================

function showFileName() {

    resetWorkflow();

    const files =
        getFiles();


    setText(

        "fileName",

        files.length

            ? `${files.length} FASTA file${files.length > 1 ? "s" : ""} selected`

            : "No files selected"

    );


    const selected =
        $("selectedFiles");


    if (!selected) return;


    selected.innerHTML =
        files.map(
            (file, index) => `

                <div class="file-chip">

                    ${index + 1}.
                    ${escapeHTML(file.name)}

                </div>

            `
        ).join("");

}


// ============================================================
// JSON
// ============================================================

async function parseJSONResponse(
    response
) {

    const data =
        await response.json();


    if (
        !response.ok
        ||
        !data.success
    ) {

        throw new Error(

            data.message

            ||

            (
                Array.isArray(
                    data.errors
                )

                ?

                data.errors.join("; ")

                :

                "Request failed."
            )

        );

    }

    return data;

}


// ============================================================
// VALIDATION
// ============================================================

async function validateFiles() {

    const files =
        getFiles();


    const status =
        $("validationStatus");


    if (!files.length) {

        showError(

            status,

            "Validation Failed",

            "Please select at least one FASTA file."

        );

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
                    body:
                        createFormData()
                }
            );


        const data =
            await parseJSONResponse(
                response
            );


        const stats =
            data.statistics || {};


        setText(
            "sequenceCount",
            stats.number_of_sequences
        );


        setText(
            "sequenceType",
            stats.sequence_type
        );


        setText(
            "totalLength",
            stats.total_length
        );


        setText(
            "fileCount",
            stats.number_of_files
        );


        setText(
            "summarySequences",
            stats.number_of_sequences
        );


        showSuccess(

            status,

            "FASTA Validation Successful ✓",

            `${stats.number_of_sequences} valid sequences detected.`

        );


        validated = true;

        $("alignBtn").disabled = false;


    }

    catch (error) {

        showError(

            status,

            "Validation Failed",

            error.message

        );

    }

}


// ============================================================
// ALIGNMENT
// ============================================================

async function runAlignment() {

    if (!validated) {

        showError(

            $("alignmentResult"),

            "Alignment Not Ready",

            "Validate the FASTA files first."

        );

        return;

    }


    const preview =
        $("alignmentPreview");


    $("alignBtn").disabled = true;


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
                    body:
                        createFormData()
                }
            );


        const data =
            await parseJSONResponse(
                response
            );


        window.phyloAlignment =
            data.alignment;


        alignmentReady = true;


        renderAlignment(
            data.alignment
        );


        setText(

            "alignmentMeta",

            `${data.number_of_sequences} sequences | ${data.alignment_length} positions | ${data.execution_time} seconds`

        );


        setText(

            "summaryAlignment",

            `${data.alignment_length} positions`

        );


        $("alignmentResult").innerHTML = `

            <div class="result-card">

                <h3>
                    Multiple Sequence Alignment Completed ✓
                </h3>

                <p>

                    <strong>
                        Sequences:
                    </strong>

                    ${escapeHTML(
                        data.number_of_sequences
                    )}

                </p>

                <p>

                    <strong>
                        Alignment length:
                    </strong>

                    ${escapeHTML(
                        data.alignment_length
                    )}

                </p>

                <p>

                    <strong>
                        Execution time:
                    </strong>

                    ${escapeHTML(
                        data.execution_time
                    )}

                    seconds

                </p>

            </div>

        `;


        if (
            data.number_of_sequences >= 3
        ) {

            $("njBtn").disabled = false;

            $("mlBtn").disabled = false;

        }


    }

    catch (error) {

        showError(

            preview,

            "Alignment Error",

            error.message

        );

    }

    finally {

        $("alignBtn").disabled =
            !validated;

    }

}


// ============================================================
// ALIGNMENT DISPLAY
// ============================================================

function getBaseClass(
    character,
    conserved
) {

    let cls = {

        A: "base-A",

        T: "base-T",

        G: "base-G",

        C: "base-C",

        U: "base-U",

        "-": "base-gap"

    }[character]

    ||

    "base-other";


    if (conserved) {

        cls +=
            " conserved";

    }


    return cls;

}


function renderAlignment(rows) {

    const preview =
        $("alignmentPreview");


    if (
        !rows?.length
    ) {

        preview.innerHTML =
            "<p>No alignment available.</p>";

        return;

    }


    const sequences =
        rows.map(
            r =>
                String(
                    r.sequence || ""
                ).toUpperCase()
        );


    const length =
        Math.max(
            ...sequences.map(
                s => s.length
            )
        );


    let consensus = "";


    for (
        let p = 0;
        p < length;
        p++
    ) {

        const counts = {};


        rows.forEach(
            row => {

                const c =
                    (
                        row.sequence[p]
                        ||
                        "-"
                    ).toUpperCase();


                if (c !== "-") {

                    counts[c] =
                        (
                            counts[c]
                            ||
                            0
                        ) + 1;

                }

            }
        );


        const best =
            Object.entries(
                counts
            ).sort(
                (a, b) =>
                    b[1] - a[1]
            )[0];


        consensus +=
            best
                ? best[0]
                : "-";

    }


    const consensusHTML =
        [...consensus]
        .map(
            (c, i) => `

                <span
                    class="consensus-char"
                    title="Position ${i + 1}"
                >

                    ${escapeHTML(c)}

                </span>

            `
        )
        .join("");


    const rowsHTML =
        rows.map(
            row => {

                const seq =
                    String(
                        row.sequence || ""
                    ).toUpperCase();


                const chars = [];


                for (
                    let i = 0;
                    i < length;
                    i++
                ) {

                    const c =
                        seq[i]
                        ||
                        "-";


                    chars.push(`

                        <span
                            class="
                                base
                                ${getBaseClass(
                                    c,
                                    c !== "-"
                                    &&
                                    c === consensus[i]
                                )}
                            "
                            title="
                                Position ${i + 1}
                            "
                        >

                            ${escapeHTML(c)}

                        </span>

                    `);

                }


                return `

                    <div class="alignment-row">

                        <div
                            class="sequence-name"
                            title="${escapeHTML(
                                row.name
                            )}"
                        >

                            ${escapeHTML(
                                row.name
                            )}

                        </div>


                        <div class="sequence-track">

                            ${chars.join("")}

                        </div>

                    </div>

                `;

            }
        )
        .join("");


    preview.innerHTML = `

        <div class="alignment-container">

            <div
                class="
                    alignment-row
                    consensus-row
                "
            >

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
// NJ TREE
// ============================================================

async function buildNJTree() {

    if (!alignmentReady)
        return;


    $("njBtn").disabled = true;


    $("njStatus").innerHTML = `

        <div class="loading compact">

            Building Neighbor Joining tree...

        </div>

    `;


    try {

        const response =
            await fetch(
                "/build-nj",
                {
                    method: "POST",
                    body:
                        createFormData()
                }
            );


        const data =
            await parseJSONResponse(
                response
            );


        phyloTreeData =
            {
                ...(phyloTreeData || {}),
                ...data
            };


        njReady = true;


        $("njNewick").textContent =
            data.nj_newick;


        showSuccess(

            $("njStatus"),

            "NJ Tree Ready ✓",

            `${data.execution_time} seconds.`

        );


        updateHybridButton();


    }

    catch (error) {

        showError(

            $("njStatus"),

            "NJ Tree Error",

            error.message

        );

    }

    finally {

        $("njBtn").disabled =
            false;

    }

}


// ============================================================
// ML TREE
// ============================================================

async function buildMLTree() {

    if (!alignmentReady)
        return;


    $("mlBtn").disabled = true;


    $("mlStatus").innerHTML = `

        <div class="loading compact">

            Building Maximum Likelihood tree...

        </div>

    `;


    try {

        const response =
            await fetch(
                "/build-ml",
                {
                    method: "POST",
                    body:
                        createFormData()
                }
            );


        const data =
            await parseJSONResponse(
                response
            );


        phyloTreeData =
            {
                ...(phyloTreeData || {}),
                ...data
            };


        mlReady = true;


        $("mlNewick").textContent =
            data.ml_newick;


        showSuccess(

            $("mlStatus"),

            "ML Tree Ready ✓",

            `${data.execution_time} seconds.`

        );


        updateHybridButton();


    }

    catch (error) {

        showError(

            $("mlStatus"),

            "ML Tree Error",

            error.message

        );

    }

    finally {

        $("mlBtn").disabled =
            false;

    }

}


function updateHybridButton() {

    $("hybridBtn").disabled =
        !(
            njReady
            &&
            mlReady
        );

}


// ============================================================
// HYBRID TREE
// ============================================================

async function buildHybridTree() {

    if (
        !(
            njReady
            &&
            mlReady
        )
    )
        return;


    $("hybridBtn").disabled =
        true;


    $("hybridStatus").innerHTML = `

        <div class="loading">

            Generating hybrid tree:

            NJ starting topology

            →

            ML refinement...

        </div>

    `;


    try {

        const response =
            await fetch(
                "/build-hybrid",
                {
                    method: "POST",
                    body:
                        createFormData()
                }
            );


        const data =
            await parseJSONResponse(
                response
            );


        phyloTreeData =
            data;


        hybridReady = true;


        $("hybridNewick").textContent =
            data.final_newick;


        showSuccess(

            $("hybridStatus"),

            "Hybrid Tree Generated ✓",

            "NJ starting topology was refined using Maximum Likelihood."

        );


        $("bootstrapBtn").disabled =
            false;


        setText(
            "summaryTree",
            "Hybrid Tree Ready"
        );


        renderInteractiveTree(
            data.final_newick,
            null
        );


        $("visualization-section")
            .scrollIntoView({
                behavior: "smooth"
            });


    }

    catch (error) {

        showError(

            $("hybridStatus"),

            "Hybrid Tree Error",

            error.message

        );

    }

    finally {

        $("hybridBtn").disabled =
            false;

    }

}


// ============================================================
// ADAPTIVE BOOTSTRAP
// ============================================================

async function runAdaptiveBootstrap() {

    if (!hybridReady) {

        showError(

            $("bootstrapResult"),

            "Bootstrap Not Ready",

            "Generate the hybrid tree first."

        );

        return;

    }


    const initial =
        Math.round(
            clamp(
                getNumberValue(
                    "initialReplicates",
                    100
                ),
                10,
                1000
            )
        );


    const batch =
        Math.round(
            clamp(
                getNumberValue(
                    "batchSize",
                    100
                ),
                10,
                500
            )
        );


    const maximum =
        Math.round(
            clamp(
                getNumberValue(
                    "maxReplicates",
                    1000
                ),
                initial,
                2000
            )
        );


    const threshold =
        clamp(
            getNumberValue(
                "stabilityThreshold",
                2
            ),
            0.1,
            20
        );


    $("bootstrapBtn").disabled =
        true;


    $("bootstrapRound1").innerHTML =
        "";

    $("bootstrapRound2").innerHTML =
        "";

    $("bootstrapResult").innerHTML =
        "";


    $("bootstrapStatus").innerHTML = `

        <div class="loading">

            Round 1:

            calculating initial branch confidence

            using ${initial} replicates...

        </div>

    `;


    try {

        // ----------------------------------------------------
        // ROUND 1
        // ----------------------------------------------------

        const fd =
            createFormData();


        fd.append(
            "initial_replicates",
            initial
        );


        fd.append(
            "batch_size",
            batch
        );


        fd.append(
            "max_replicates",
            maximum
        );


        fd.append(
            "stability_threshold",
            threshold
        );


        const response1 =
            await fetch(
                "/bootstrap-round1",
                {
                    method: "POST",
                    body: fd
                }
            );


        const round1 =
            await parseJSONResponse(
                response1
            );


        bootstrapToken =
            round1.token;


        bootstrapData =
            round1;


        // IMPORTANT:
        // Show Round 1 before Round 2
        renderBootstrapRound1(
            round1
        );


        // Visualize initial confidence
        renderInteractiveTree(

            round1.hybrid_newick,

            round1.branches

        );


        // ----------------------------------------------------
        // STOP IF ALREADY STABLE
        // ----------------------------------------------------

        if (
            !round1.needs_round2
        ) {

            finishBootstrap(

                round1,

                "Round 1 was sufficient; no additional replicates were required."

            );

            return;

        }


        // ----------------------------------------------------
        // ROUND 2
        // ----------------------------------------------------

        $("bootstrapStatus").innerHTML = `

            <div class="loading">

                Round 1 completed.

                Additional sampling is required.

                Running Round 2 automatically...

            </div>

        `;


        const fd2 =
            new FormData();


        fd2.append(
            "token",
            bootstrapToken
        );


        const response2 =
            await fetch(

                "/bootstrap-round2",

                {
                    method: "POST",
                    body: fd2
                }

            );


        const round2 =
            await parseJSONResponse(
                response2
            );


        bootstrapData =
            round2;


        renderBootstrapRound2(
            round2
        );


        renderInteractiveTree(

            round2.hybrid_newick,

            round2.branches

        );


        finishBootstrap(

            round2,

            "Adaptive bootstrap completed."

        );


    }

    catch (error) {

        showError(

            $("bootstrapResult"),

            "Adaptive Bootstrap Error",

            error.message

        );

        $("bootstrapStatus")
            .innerHTML = "";

    }

    finally {

        $("bootstrapBtn").disabled =
            false;

    }

}


// ============================================================
// BOOTSTRAP TABLE
// ============================================================

function renderBranchTable(
    branches
) {

    if (
        !Array.isArray(branches)
        ||
        !branches.length
    ) {

        return `

            <p class="muted">

                No internal branches
                available.

            </p>

        `;

    }


    return `

        <div class="branch-table-wrap">

            <table class="branch-table">

                <thead>

                    <tr>

                        <th>
                            Branch
                        </th>

                        <th>
                            Confidence
                        </th>

                        <th>
                            Class
                        </th>

                    </tr>

                </thead>


                <tbody>

                    ${branches.map(
                        b => `

                            <tr>

                                <td>

                                    ${escapeHTML(
                                        b.branch
                                    )}

                                </td>


                                <td>

                                    <strong>

                                        ${escapeHTML(
                                            b.support
                                        )}%

                                    </strong>

                                </td>


                                <td>

                                    <span
                                        class="
                                            confidence-badge
                                            ${b.confidence.toLowerCase()}
                                        "
                                    >

                                        ${escapeHTML(
                                            b.confidence
                                        )}

                                    </span>

                                </td>

                            </tr>

                        `
                    ).join("")}

                </tbody>

            </table>

        </div>

    `;

}


// ============================================================
// ROUND 1 DISPLAY
// ============================================================

function renderBootstrapRound1(
    data
) {

    $("bootstrapRound1").innerHTML = `

        <div
            class="
                round-card
                initial
            "
        >

            <div class="round-header">

                <div>

                    <span class="round-pill">
                        ROUND 1
                    </span>

                    <h3>
                        Initial Bootstrap Results
                    </h3>

                </div>


                <strong>

                    ${escapeHTML(
                        data.replicates
                    )}

                    replicates

                </strong>

            </div>


            <p>

                Initial branch confidence
                is shown before the adaptive
                decision.

            </p>


            <div
                class="
                    result-grid
                    compact-grid
                "
            >

                <div class="result-item">

                    <span>
                        Average Confidence
                    </span>

                    <strong>

                        ${escapeHTML(
                            data.average_confidence
                        )}%

                    </strong>

                </div>


                <div class="result-item">

                    <span>
                        Branches
                    </span>

                    <strong>

                        ${escapeHTML(
                            data.branches.length
                        )}

                    </strong>

                </div>


                <div class="result-item">

                    <span>
                        Additional Sampling
                    </span>

                    <strong>

                        ${
                            data.needs_round2
                            ? "Required"
                            : "Not required"
                        }

                    </strong>

                </div>

            </div>


            ${renderBranchTable(
                data.branches
            )}

        </div>

    `;

}


// ============================================================
// ROUND 2 DISPLAY
// ============================================================

function renderBootstrapRound2(
    data
) {

    $("bootstrapRound2").innerHTML = `

        <div
            class="
                round-card
                final-round
            "
        >

            <div class="round-header">

                <div>

                    <span class="round-pill">
                        ROUND 2
                    </span>

                    <h3>
                        Adaptive Bootstrap Completed
                    </h3>

                </div>


                <strong>

                    ${escapeHTML(
                        data.replicates
                    )}

                    total replicates

                </strong>

            </div>


            <p>

                Additional sampling was
                performed because the initial
                analysis required more evidence.

            </p>


            ${renderBranchTable(
                data.branches
            )}

        </div>

    `;

}


// ============================================================
// FINAL RESULT
// ============================================================

function finishBootstrap(
    data,
    message
) {

    bootstrapReady = true;


    const unreliable =
        Array.isArray(
            data.unreliable_regions
        )

        ?

        data.unreliable_regions

        :

        [];


    $("bootstrapStatus").innerHTML = `

        <div class="success-box">

            <h3>
                ✓ Adaptive Bootstrap Completed
            </h3>

            <p>
                ${escapeHTML(message)}
            </p>

        </div>

    `;


    // FINAL RESULT ONLY
    $("bootstrapResult").innerHTML = `

        <div
            class="
                result-card
                final-bootstrap-result
            "
        >

            <h3>
                Final Adaptive Bootstrap Result ✓
            </h3>


            <div class="result-grid">


                <div class="result-item">

                    <span>
                        Total Replicates
                    </span>

                    <strong>

                        ${escapeHTML(
                            data.replicates
                        )}

                    </strong>

                </div>


                <div class="result-item">

                    <span>
                        Average Confidence
                    </span>

                    <strong>

                        ${escapeHTML(
                            data.average_confidence
                        )}%

                    </strong>

                </div>


                <div class="result-item">

                    <span>
                        Stability
                    </span>

                    <strong
                        class="
                            ${
                                data.stable
                                ? "success-text"
                                : ""
                            }
                        "
                    >

                        ${
                            data.stable
                            ? "Stable"
                            : "Maximum reached"
                        }

                    </strong>

                </div>


                <div class="result-item">

                    <span>
                        Unreliable Regions
                    </span>

                    <strong>

                        ${escapeHTML(
                            unreliable.length
                        )}

                    </strong>

                </div>


            </div>


            <p>

                <strong>
                    Stopping reason:
                </strong>

                ${escapeHTML(
                    data.stopping_reason
                    ||
                    "Bootstrap completed"
                )}

            </p>


            ${
                unreliable.length

                ?

                `

                <div
                    class="
                        bootstrap-warning
                    "
                >

                    <h4>
                        ⚠ Unreliable Regions
                    </h4>

                    <ul>

                        ${
                            unreliable
                            .map(
                                b => `

                                    <li>

                                        ${escapeHTML(
                                            b.branch
                                        )}

                                        —

                                        ${escapeHTML(
                                            b.support
                                        )}%

                                    </li>

                                `
                            )
                            .join("")
                        }

                    </ul>

                </div>

                `

                :

                `

                <div class="success-box">

                    <strong>

                        No unreliable regions detected.

                    </strong>

                </div>

                `
            }


        </div>

    `;


    setText(
        "summaryBootstrap",
        `${data.replicates} reps`
    );


    $("finalResults").innerHTML = `

        <div class="result-card">

            <h3>
                Project Outcome
            </h3>

            <p>

                <strong>
                    Tree:
                </strong>

                Hybrid NJ + ML

            </p>

            <p>

                <strong>
                    Bootstrap:
                </strong>

                Adaptive

            </p>

            <p>

                <strong>
                    Final confidence:
                </strong>

                ${escapeHTML(
                    data.average_confidence
                )}%

            </p>

            <p>

                <strong>
                    Unreliable regions:
                </strong>

                ${escapeHTML(
                    unreliable.length
                )}

            </p>

        </div>

    `;

}


// ============================================================
// NEWICK PARSER
// ============================================================

function tokenizeNewick(
    text
) {

    return text

        .replace(
            /;/g,
            ""
        )

        .match(
            /\(|\)|,|:[^(),;]+|[^(),;:]+/g
        )

        || [];

}


function parseNewick(
    text
) {

    const tokens =
        tokenizeNewick(
            text.trim()
        );


    let pos = 0;


    function parseNode() {

        if (
            tokens[pos] === "("
        ) {

            pos++;


            const children = [
                parseNode()
            ];


            while (
                tokens[pos] === ","
            ) {

                pos++;

                children.push(
                    parseNode()
                );

            }


            if (
                tokens[pos] === ")"
            ) {

                pos++;

            }


            let label = "";

            let length = 0;


            if (
                tokens[pos]
                &&
                ![
                    "(",
                    ")",
                    ","
                ].includes(
                    tokens[pos]
                )
            ) {

                const token =
                    tokens[pos++];


                if (
                    token.startsWith(":")
                ) {

                    length =
                        Number(
                            token.slice(1)
                        )
                        ||
                        0;

                }

                else {

                    label = token;


                    if (
                        tokens[pos]
                        ?.startsWith(":")
                    ) {

                        length =
                            Number(
                                tokens[pos++]
                                    .slice(1)
                            )
                            ||
                            0;

                    }

                }

            }


            return {

                name: label,

                length: length,

                children: children

            };

        }


        let name =
            tokens[pos++]
            ||
            "Unknown";


        let length = 0;


        if (
            tokens[pos]
            ?.startsWith(":")
        ) {

            length =
                Number(
                    tokens[pos++]
                        .slice(1)
                )
                ||
                0;

        }


        return {

            name,

            length,

            children: []

        };

    }


    return parseNode();

}


function normalizeName(
    name
) {

    return String(
        name || ""
    )
    .trim()
    .replace(
        /^['"]|['"]$/g,
        ""
    );

}


function descendantNames(
    node
) {

    if (
        !node.children?.length
    ) {

        return [
            normalizeName(
                node.name
            )
        ];

    }


    return node.children
        .flatMap(
            descendantNames
        )
        .filter(Boolean);

}


function confidenceForNode(
    node,
    branches
) {

    if (
        !branches?.length
        ||
        !node.children?.length
    ) {

        return null;

    }


    const names =
        descendantNames(node)
            .sort();


    return branches.find(
        b =>

            Array.isArray(
                b.members
            )

            &&

            b.members
                .slice()
                .sort()
                .join("|")

            ===

            names.join("|")

    )
    ||
    null;

}


// ============================================================
// INTERACTIVE TREE
// ============================================================

function confidenceClass(
    value
) {

    if (value == null)
        return "neutral";


    if (value >= 80)
        return "high";


    if (value >= 50)
        return "moderate";


    return "low";

}


function renderInteractiveTree(
    newick,
    branches
) {

    const canvas =
        $("treeCanvas");


    if (
        !canvas
        ||
        !newick
    )
        return;


    let root;


    try {

        root =
            parseNewick(
                newick
            );

    }

    catch {

        canvas.innerHTML = `

            <div class="empty-viz">

                Could not parse
                the Newick tree.

            </div>

        `;

        return;

    }


    const leaves = [];


    function collect(
        node
    ) {

        if (
            !node.children?.length
        ) {

            leaves.push(
                node
            );

        }

        else {

            node.children.forEach(
                collect
            );

        }

    }


    collect(root);


    const rowHeight = 48;

    const width = 1050;

    const height =
        Math.max(
            380,
            leaves.length
            *
            rowHeight
            +
            80
        );


    let leafIndex = 0;


    function assignY(
        node
    ) {

        if (
            !node.children?.length
        ) {

            node.y =
                40
                +
                leafIndex++
                *
                rowHeight;


            return node.y;

        }


        node.children.forEach(
            assignY
        );


        node.y =
            node.children.reduce(
                (
                    sum,
                    child
                ) =>
                    sum + child.y,
                0
            )
            /
            node.children.length;


        return node.y;

    }


    assignY(root);


    function maxDepth(
        node,
        depth = 0
    ) {

        if (
            !node.children?.length
        )
            return depth;


        return Math.max(

            ...node.children.map(

                child =>
                    maxDepth(
                        child,
                        depth + 1
                    )

            )

        );

    }


    const depth =
        maxDepth(root);


    const xStep =
        Math.max(
            110,
            Math.min(
                190,
                (
                    width - 220
                )
                /
                Math.max(
                    depth,
                    1
                )
            )
        );


    function assignX(
        node,
        x = 50
    ) {

        node.x = x;


        if (
            node.children?.length
        ) {

            node.children.forEach(
                child =>
                    assignX(
                        child,
                        x + xStep
                    )
            );

        }

    }


    assignX(root);


    let branchesSVG = "";

    let labelsSVG = "";


    function draw(
        node
    ) {

        if (
            !node.children?.length
        ) {

            labelsSVG += `

                <text
                    x="${node.x + 10}"
                    y="${node.y + 5}"
                    class="leaf-label"
                >

                    ${escapeHTML(
                        normalizeName(
                            node.name
                        )
                    )}

                </text>

            `;

            return;

        }


        const branch =
            confidenceForNode(
                node,
                branches
            );


        const cls =
            confidenceClass(
                branch?.support
            );


        const colors = {

            high: "#238b70",

            moderate: "#d18b1f",

            low: "#d44a4a",

            neutral: "#7a8798"

        };


        const color =
            colors[cls];


        const childrenY =
            node.children.map(
                child =>
                    child.y
            );


        branchesSVG += `

            <line

                x1="${node.x}"

                y1="${Math.min(
                    ...childrenY
                )}"

                x2="${node.x}"

                y2="${Math.max(
                    ...childrenY
                )}"

                class="
                    tree-branch
                    ${cls}
                "

                style="
                    stroke:${color}
                "

            ></line>

        `;


        node.children.forEach(
            child => {

                branchesSVG += `

                    <line

                        x1="${node.x}"

                        y1="${node.y}"

                        x2="${child.x}"

                        y2="${child.y}"

                        class="
                            tree-branch
                            ${cls}
                        "

                        style="
                            stroke:${color}
                        "

                    ></line>

                `;


                draw(child);

            }
        );


        node._id =
            `n${
                Math.random()
                .toString(36)
                .slice(2)
            }`;


        const supportText =
            branch
                ? `${branch.support}%`
                : "";


        branchesSVG += `

            <circle

                cx="${node.x}"

                cy="${node.y}"

                r="8"

                class="
                    tree-node
                    ${cls}
                "

                data-node-id="
                    ${node._id}
                "

                tabindex="0"

            ></circle>


            ${
                supportText

                ?

                `

                <text

                    x="${node.x - 8}"

                    y="${node.y - 14}"

                    class="support-label"

                >

                    ${escapeHTML(
                        supportText
                    )}

                </text>

                `

                :

                ""
            }

        `;

    }


    draw(root);


    canvas.innerHTML = `

        <svg

            id="treeSvg"

            viewBox="
                0
                0
                ${width}
                ${height}
            "

            role="img"

            aria-label="
                Interactive phylogenetic tree
            "

        >

            <g
                id="treeGroup"
                transform="
                    scale(${treeZoom})
                "
            >

                ${branchesSVG}

                ${labelsSVG}

            </g>

        </svg>

    `;


    function attach(
        node
    ) {

        if (
            node.children?.length
        ) {

            const branch =
                confidenceForNode(
                    node,
                    branches
                );


            const nodeEl =
                canvas.querySelector(
                    `[data-node-id="${node._id}"]`
                );


            if (nodeEl) {

                nodeEl.addEventListener(
                    "click",
                    () =>
                        showBranchDetails(
                            node,
                            branch
                        )
                );

            }


            node.children.forEach(
                attach
            );

        }

    }


    attach(root);


    setText(

        "vizHint",

        branches?.length

        ?

        "Click an internal node to inspect bootstrap confidence."

        :

        "Hybrid tree displayed. Run bootstrap to add confidence."

    );

}


// ============================================================
// BRANCH DETAILS
// ============================================================

function showBranchDetails(
    node,
    branch
) {

    const details =
        $("branchDetails");


    const members =
        descendantNames(
            node
        );


    if (!branch) {

        details.innerHTML = `

            <strong>
                Internal branch
            </strong>

            <p>

                Confidence has not
                been calculated yet.

            </p>

            <p>

                <strong>
                    Descendants:
                </strong>

                ${escapeHTML(
                    members.join(", ")
                )}

            </p>

        `;

        return;

    }


    details.innerHTML = `

        <div class="detail-title">

            Selected Branch

        </div>


        <div class="detail-value">

            ${escapeHTML(
                branch.support
            )}%

            —

            ${escapeHTML(
                branch.confidence
            )}

            confidence

        </div>


        <p>

            <strong>
                Descendants:
            </strong>

            ${escapeHTML(
                members.join(", ")
            )}

        </p>

    `;

}


// ============================================================
// ZOOM
// ============================================================

function applyTreeZoom() {

    const group =
        $("treeGroup");


    if (group) {

        group.setAttribute(

            "transform",

            `scale(${treeZoom})`

        );

    }

}


function setupZoom() {

    $("zoomInBtn")
        ?.addEventListener(
            "click",
            () => {

                treeZoom =
                    clamp(
                        treeZoom + 0.1,
                        0.6,
                        1.8
                    );

                applyTreeZoom();

            }
        );


    $("zoomOutBtn")
        ?.addEventListener(
            "click",
            () => {

                treeZoom =
                    clamp(
                        treeZoom - 0.1,
                        0.6,
                        1.8
                    );

                applyTreeZoom();

            }
        );


    $("zoomResetBtn")
        ?.addEventListener(
            "click",
            () => {

                treeZoom = 1;

                applyTreeZoom();

            }
        );

}


// ============================================================
// STARTUP
// ============================================================

document.addEventListener(
    "DOMContentLoaded",
    () => {

        $("fastaFile")
            ?.addEventListener(
                "change",
                showFileName
            );


        $("validateBtn")
            ?.addEventListener(
                "click",
                validateFiles
            );


        $("alignBtn")
            ?.addEventListener(
                "click",
                runAlignment
            );


        $("njBtn")
            ?.addEventListener(
                "click",
                buildNJTree
            );


        $("mlBtn")
            ?.addEventListener(
                "click",
                buildMLTree
            );


        $("hybridBtn")
            ?.addEventListener(
                "click",
                buildHybridTree
            );


        $("bootstrapBtn")
            ?.addEventListener(
                "click",
                runAdaptiveBootstrap
            );


        setupZoom();

    }
);