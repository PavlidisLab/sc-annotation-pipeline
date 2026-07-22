pipeline {
    agent any

    stages {
        stage('Test: mouse') {
            steps {
                sh '''
                    set -euo pipefail
                    nextflow run main.nf -profile test_mouse,singularity --upload_gemma false -work-dir work/ci_test_mouse -process.executor slurm -resume
                '''
            }
        }
        stage('Test: human') {
            steps {
                sh '''
                    set -euo pipefail
                    nextflow run main.nf -profile test_human,singularity --upload_gemma false -work-dir work/ci_test_human -process.executor slurm -resume
                '''
            }
        }
        stage('Test: mouse (per-sample)') {
            steps {
                sh '''
                    set -euo pipefail
                    nextflow run main.nf -profile test_mouse_persample,singularity --upload_gemma false -work-dir work/ci_test_mouse_persample -process.executor slurm -resume
                '''
            }
        }
        stage('Test: mouse (local study_paths)') {
            steps {
                sh '''
                    set -euo pipefail
                    nextflow run main.nf -profile test_mouse_local,singularity --upload_gemma false -work-dir work/ci_test_mouse_local -process.executor slurm -resume
                '''
            }
        }
        stage('Test: human (local study_paths)') {
            steps {
                sh '''
                    set -euo pipefail
                    nextflow run main.nf -profile test_human_local,singularity --upload_gemma false -work-dir work/ci_test_human_local -process.executor slurm -resume
                '''
            }
        }
        stage('Test: mouse (Gemma upload disabled)') {
            steps {
                sh '''
                    set -euo pipefail
                    nextflow run main.nf -profile test_mouse_upload_off,singularity --upload_gemma false -work-dir work/ci_test_mouse_upload_off -process.executor slurm -resume
                '''
            }
        }
    }

    post {
        success {
            echo 'All test profiles passed - safe to tag a release.'
        }
        failure {
            echo 'One or more test profiles failed - do not tag a release.'
        }
    }
}
