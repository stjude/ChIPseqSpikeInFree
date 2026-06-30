param(
  [string]$RepoRoot = (Get-Location).Path
)

$ErrorActionPreference = "Stop"

$rscriptCmd = Get-Command Rscript -ErrorAction SilentlyContinue
if ($rscriptCmd) {
  $rscript = $rscriptCmd.Source
} else {
  $fallbacks = @(
    "C:\Users\hjin\AppData\Local\Programs\R\R-4.3.3\bin\x64\Rscript.exe",
    "C:\Users\hjin\AppData\Local\Programs\R\R-4.3.3\bin\Rscript.exe"
  )
  $rscript = $fallbacks | Where-Object { Test-Path $_ } | Select-Object -First 1
}

if (-not $rscript) {
  throw "Rscript was not found on PATH or in the known Windows fallback locations."
}

Push-Location $RepoRoot
try {
  & $rscript -e "testthat::test_dir('tests/testthat')"
} finally {
  Pop-Location
}