graph TD
    %% Inputs
    A((study_names))
    B((study_paths))
    C((params.* and files))

    %% Processes
    P1[INPUT_CHECK]
    P2[PREPARE_REFERENCE]
    P3[PROCESS_QUERIES]
    P4[CLASSIFY_CELLTYPES]
    P5[QC_REPORTING]
    %% Uncomment if GEMMA_UPLOAD is enabled
    %% P6[GEMMA_UPLOAD]

    %% Outputs
    O1((celltypes))
    O2((masks))
    O3((multiqc))
    O4((versions))

    %% Flow
    A --> P1
    B --> P1
    C --> P2
    P1 -->|studies| P3
    P2 -->|model_path| P3
    P2 -->|ref_paths| P4
    P3 -->|processed| P4
    P3 -->|raw| P5
    P1 -->|versions| O4
    P2 -->|versions| O4
    P3 -->|versions| O4
    P4 -->|celltypes| P5
    P4 -->|celltypes| O1
    P4 -->|versions| O4
    P5 -->|masks| O2
    P5 -->|multiqc| O3
    P5 -->|versions| O4
    %% Uncomment if GEMMA_UPLOAD is enabled
    %% P4 --> P6
    %% P5 --> P6
    %% P6 --> O5((gemma_upload_message))
