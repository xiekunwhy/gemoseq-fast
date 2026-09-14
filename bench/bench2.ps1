param(
    [Parameter(Mandatory=$true)][string]$Bam,
    [Parameter(Mandatory=$true)][string]$Prefix,
    [string]$Jar = "GeMoSeq-1.2.3-fast.jar",
    [string]$Genome = "genome.fa",
    [string]$Strand = "FR_UNSTRANDED",
    [string]$Xmx = "16g",
    [int]$Threads = 6,
    [string]$Extra = "",
    [int]$TimeoutSec = 0,
    [string]$WorkDir = "",
    [string]$JavaExe = "C:\Program Files\Java\jdk-11.0.11\bin\java.exe"
)
if (-not (Test-Path $JavaExe)) { $JavaExe = "java" }
$dir = if ($WorkDir -ne "") { $WorkDir } else { Split-Path -Parent $MyInvocation.MyCommand.Path }
Set-Location $dir
$outdir = "${Prefix}_out"
New-Item -ItemType Directory -Force -Path $outdir | Out-Null
$stdout = Join-Path $dir "${Prefix}.stdout.log"
$stderr = Join-Path $dir "${Prefix}.stderr.log"
$argList = @("-Xmx$Xmx",'-jar',$Jar,'gemoseq',"g=$Genome","m=$Bam","s=$Strand","threads=$Threads","outdir=$outdir")
if ($Extra -ne "") { $argList += $Extra -split '\s+' }
$argStr = ($argList | ForEach-Object { if ($_ -match '\s') { '"' + $_ + '"' } else { $_ } }) -join ' '
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $JavaExe
$psi.Arguments = $argStr
$psi.WorkingDirectory = $dir
$psi.UseShellExecute = $false
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$p = New-Object System.Diagnostics.Process
$p.StartInfo = $psi
$t0 = Get-Date
[void]$p.Start()
$outTask = $p.StandardOutput.ReadToEndAsync()
$errTask = $p.StandardError.ReadToEndAsync()
$peak = 0
$killed = $false
while (-not $p.HasExited) {
    Start-Sleep -Milliseconds 500
    try { $p.Refresh(); if ($p.PeakWorkingSet64 -gt $peak) { $peak = $p.PeakWorkingSet64 } } catch {}
    if ($TimeoutSec -gt 0 -and ((Get-Date) - $t0).TotalSeconds -gt $TimeoutSec) {
        try { $p.Kill() } catch {}
        $killed = $true
        break
    }
}
try { $p.WaitForExit() } catch {}
$wall = [int](((Get-Date) - $t0).TotalSeconds)
$code = -999
try { $code = $p.ExitCode } catch {}
try { [System.IO.File]::WriteAllText($stdout, $outTask.Result) } catch {}
try { [System.IO.File]::WriteAllText($stderr, $errTask.Result) } catch {}
$peakMB = [int]($peak / 1MB)
$res = "bam=$Bam jar=$Jar genome=$Genome strand=$Strand extra=$Extra exit=$code peakMB=$peakMB wall=${wall}s killed=$killed"
Write-Output $res
Set-Content -Path (Join-Path $dir "${Prefix}.bench.txt") -Value $res
