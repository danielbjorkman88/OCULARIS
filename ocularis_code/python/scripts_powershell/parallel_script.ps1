$patientIDs = "P1", "P2", "P3"  # List of patient IDs
$margins = 0.1, 0.2, 0.3  # List of margin values

# Iterate over patient IDs
foreach ($patientID in $patientIDs) {
    # Iterate over margins
    foreach ($margin in $margins) {
        # Define the command to run the Python script with patientID and margin parameters
        $command = "C:\GitHub\Optistps\python\multi_process_test.py $patientID $margin"

        # Start a job for each combination of patientID and margin
        $job = Start-Job -ScriptBlock { param($command) Invoke-Expression $command } -ArgumentList $command

        # Wait for the job to finish
        Wait-Job $job | Out-Null

        # Retrieve and display the job result
        $result = Receive-Job $job
        Write-Output $result

        # Remove the job
        Remove-Job $job
    }
}