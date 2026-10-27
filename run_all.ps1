<#
.SYNOPSIS
  Runs the experiments of the paper with unrefinement_example and compares
  the final number of MPBES functions with the logs in runs\<run>\.

.EXAMPLE
  .\run_all.ps1 -Exe C:\gismo\build\bin\Release\unrefinement_example.exe -CasadiBin C:\casadi\bin
  .\run_all.ps1 -Exe ... -Only yeti_eg006_ef006,rotor_lambda0

  The arguments of every run are those of the first line of
  runs\<run>\cmd_output.txt. The output of a run is written to
  <Out>\<run>\ (cmd_output.txt, out.log, out_final.xml).
#>
param(
    [Parameter(Mandatory = $true)][string]$Exe,   # unrefinement_example.exe
    [string]$CasadiBin = "",                      # folder with casadi.dll, if not on PATH
    [string]$Out = "output",                      # relative to this script
    [string[]]$Only = @()                         # run names; default: all
)
$ErrorActionPreference = "Stop"
$Root = $PSScriptRoot
$Exe = (Resolve-Path $Exe).Path
if ($CasadiBin) { $env:PATH = (Resolve-Path $CasadiBin).Path + ";" + $env:PATH }

# run = @(input geometry in data\geometries, arguments)
$Runs = [ordered]@{
    "yeti_eg001_ef001"     = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 0.01 -e 0.01')
    "yeti_eg003_ef003"     = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 0.03 -e 0.03')
    "yeti_eg006_ef006"     = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 0.06 -e 0.06')
    "yeti_nof_ef010"       = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 1.0 -e 0.01')
    "yeti_nof_ef030"       = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 1.0 -e 0.03')
    "yeti_nof_ef060"       = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 1.0 -e 0.06')
    "yeti_F_ef010"         = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 1.0 -e 0.01 -F 0:u0,0:u1,1:u1,2:u0,2:u1,3:u1,4:u0,4:u1,5:v0,5:v1,6:v0,6:v1,7:v0,7:v1,8:v0,8:v1,9:v1,10:u0,10:u1,11:u0,12:u0,12:u1,13:u0,14:u0,14:u1,15:v1')
    "yeti_F_ef030"         = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 1.0 -e 0.03 -F 0:u0,0:u1,1:u1,2:u0,2:u1,3:u1,4:u0,4:u1,5:v0,5:v1,6:v0,6:v1,7:v0,7:v1,8:v0,8:v1,9:v1,10:u0,10:u1,11:u0,12:u0,12:u1,13:u0,14:u0,14:u1,15:v1')
    "yeti_F_ef060"         = @('yeti_mp2_thb_approximation_fine_L3_topology.xml', '-t 3 --noplot -g 1.0 -e 0.06 -F 0:u0,0:u1,1:u1,2:u0,2:u1,3:u1,4:u0,4:u1,5:v0,5:v1,6:v0,6:v1,7:v0,7:v1,8:v0,8:v1,9:v1,10:u0,10:u1,11:u0,12:u0,12:u1,13:u0,14:u0,14:u1,15:v1')
    "mask_folded_ls_only"  = @('mask_L3_PatchOnly689_noR.xml', '-t 1 --noplot')
    "mask_folded_lo_only"  = @('mask_L3_PatchOnly689_noR.xml', '-t 2 --noplot')
    "mask_folded_nlo_full" = @('mask_L3_PatchOnly689_noR.xml', '-t 3 --noplot')
    "joystick_global"      = @('joystick_approximation_fine_L3_NLO.xml', '-t 3 --noplot -e 0.020 -F 0:u0,0:u1,1:u0,1:u1,2:u0,2:u1,3:v0,3:v1,4:v0,8:v0,8:v1,9:u1,10:u0,11:v0,11:v1,12:v1,13:u0,13:u1,14:u0,14:u1,14:v0,15:v0,15:v1,16:v0,16:v1,17:u0,17:u1,17:v0,18:u0,18:u1,19:v1,20:v0,20:v1,21:u0,22:u1,23:v0,23:v1,27:v0,28:v0,28:v1,29:u0,29:u1')
    "joystick_lambda0"     = @('joystick_approximation_fine_L3_NLO.xml', '-t 3 --noplot -e 0.020 -l 0 -F 0:u0,0:u1,1:u0,1:u1,2:u0,2:u1,3:v0,3:v1,4:v0,8:v0,8:v1,9:u1,10:u0,11:v0,11:v1,12:v1,13:u0,13:u1,14:u0,14:u1,14:v0,15:v0,15:v1,16:v0,16:v1,17:u0,17:u1,17:v0,18:u0,18:u1,19:v1,20:v0,20:v1,21:u0,22:u1,23:v0,23:v1,27:v0,28:v0,28:v1,29:u0,29:u1')
    "joystick_lambda1"     = @('joystick_approximation_fine_L3_NLO.xml', '-t 3 --noplot -e 0.020 -l 1 -F 0:u0,0:u1,1:u0,1:u1,2:u0,2:u1,3:v0,3:v1,4:v0,8:v0,8:v1,9:u1,10:u0,11:v0,11:v1,12:v1,13:u0,13:u1,14:u0,14:u1,14:v0,15:v0,15:v1,16:v0,16:v1,17:u0,17:u1,17:v0,18:u0,18:u1,19:v1,20:v0,20:v1,21:u0,22:u1,23:v0,23:v1,27:v0,28:v0,28:v1,29:u0,29:u1')
    "joystick_grendas"     = @('joystick_approximation_fine_L3_NLO.xml', '-t 3 --noplot -e 0.020 -F 0:u0,0:u1,1:u0,1:u1,2:u0,2:u1,3:v0,3:v1,4:v0,8:v0,8:v1,9:u1,10:u0,11:v0,11:v1,12:v1,13:u0,13:u1,14:u0,14:u1,14:v0,15:v0,15:v1,16:v0,16:v1,17:u0,17:u1,17:v0,18:u0,18:u1,19:v1,20:v0,20:v1,21:u0,22:u1,23:v0,23:v1,27:v0,28:v0,28:v1,29:u0,29:u1')
    "joystick_lex"         = @('joystick_approximation_fine_L3_NLO.xml', '-t 3 --noplot -s lex -e 0.020 -F 0:u0,0:u1,1:u0,1:u1,2:u0,2:u1,3:v0,3:v1,4:v0,8:v0,8:v1,9:u1,10:u0,11:v0,11:v1,12:v1,13:u0,13:u1,14:u0,14:u1,14:v0,15:v0,15:v1,16:v0,16:v1,17:u0,17:u1,17:v0,18:u0,18:u1,19:v1,20:v0,20:v1,21:u0,22:u1,23:v0,23:v1,27:v0,28:v0,28:v1,29:u0,29:u1')
    "rotor_global"         = @('rotorRE.xml', '-t 3 --noplot -g 0.03 -e 0.03 -x 2')
    "rotor_lambda0"        = @('rotorRE.xml', '-t 3 --noplot -g 0.03 -e 0.03 -x 2 -l 0')
    "rotor_lambda2"        = @('rotorRE.xml', '-t 3 --noplot -g 0.03 -e 0.03 -x 2 -l 2')
    "rotor_eg0001_ef0001"  = @('rotorRE.xml', '-t 3 --noplot -g 0.001 -e 0.001 -x 2')
    "rotor_nog_ef0001"     = @('rotorRE.xml', '-t 3 --noplot -e 0.001 -x 2')
}

function Count([string]$log) {
    if (-not (Test-Path $log)) { return "-" }
    $m = [regex]::Match((Get-Content $log -Tail 1), '-> (\d+) MPBES functions')
    if ($m.Success) { return $m.Groups[1].Value } else { return "?" }
}

foreach ($name in $Runs.Keys) {
    if ($Only.Count -and -not ($Only -contains $name)) { continue }
    $geom, $arguments = $Runs[$name]
    $dir = Join-Path (Join-Path $Root $Out) $name
    New-Item -ItemType Directory -Force $dir | Out-Null
    Push-Location $dir
    try {
        $env:UNREFINEMENT_CMD = "`"$Exe`" -f `"$Root\data\geometries\$geom`" $arguments -o out"
        Set-Content -Path cmd_output.txt -Value $env:UNREFINEMENT_CMD -Encoding Ascii
        $t0 = Get-Date
        cmd /c "%UNREFINEMENT_CMD% >> cmd_output.txt 2>&1"
        $rc = $LASTEXITCODE
        $sec = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)
    } finally { Pop-Location }
    $new = Count (Join-Path $dir "out.log")
    $ref = Count (Join-Path $Root "runs\$name\out.log")
    $ok = if ($rc -eq 0 -and $new -eq $ref) { "OK" } else { "DIFFERS" }
    "{0,-24} {1,6} MPBES functions (log: {2,6})  {3,-7} {4} s" -f $name, $new, $ref, $ok, $sec
}
