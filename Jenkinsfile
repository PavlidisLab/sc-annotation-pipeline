pipeline {
    agent any

    stages {
        stage('Checkout Code') {
            steps {
                    // Checkout the repository
                checkout scm
            }
        }

        stage('Run Pipeline') {
            steps {
                    // Run the Nextflow pipeline
                sh 'nextflow run main.nf -profile conda -params-file params.mm.json'

            }
        }
    }

    post {
        success {
            echo "Pipeline ran successfully!"
        }

        failure {
            echo "Pipeline failed!"
        }
    }
}