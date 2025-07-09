pipeline {
    agent any

    stages {
        stage('Run Nextflow Pipeline') {
            steps {
                script {
                    echo 'Running sc-annotate pipeline...'

                    // Wrap in try-catch to capture errors
                    try {
                        sh '''
                            set -euo pipefail

                            nextflow run sc-annotate.nf \
                              -profile conda \
                              -params-file params.mm.json \
                              --study_names /space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/study_names_mouse.txt \
                              -process.executor slurm \
                              -resume
                        '''
                    } catch (err) {
                        echo 'Pipeline failed with the following error:'
                        echo err.getMessage()
                        error('Pipeline execution failed') // fail the Jenkins build
                    }
                }
            }
        }
    }

    post {
        success {
            echo 'Pipeline executed successfully!'
        }
        failure {
            echo 'Jenkins job failed.'
        }
    }
}
