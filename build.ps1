# SPH Particle Simulator - Windows Build Script
# Author: Lucien Vallois
# Description: Automated build script for Windows

param(
    [switch]$Debug,
    [switch]$Release,
    [switch]$NoOpenMP,
    [switch]$CUDA,
    [switch]$NoPython,
    [switch]$NoExamples,
    [switch]$NoBenchmarks,
    [switch]$Tests,
    [switch]$Verbose,
    [switch]$Clean,
    [switch]$Help
)

# Default values
$BuildType = "Release"
$BuildDir = "build"
$UseOpenMP = "ON"
$UseCUDA = "OFF"
$BuildPython = "ON"
$BuildExamples = "ON"
$BuildBenchmarks = "ON"
$BuildTests = "OFF"
$VerboseMode = "OFF"
$CleanBuild = $false

# Colors for output
$Colors = @{
    "Red" = [System.ConsoleColor]::Red
    "Green" = [System.ConsoleColor]::Green
    "Yellow" = [System.ConsoleColor]::Yellow
    "Blue" = [System.ConsoleColor]::Blue
    "White" = [System.ConsoleColor]::White
}

function Write-ColoredOutput {
    param(
        [string]$Color,
        [string]$Prefix,
        [string]$Message
    )
    $OriginalColor = $host.ui.RawUI.ForegroundColor
    $host.ui.RawUI.ForegroundColor = $Colors[$Color]
    Write-Host "[$Prefix] $Message"
    $host.ui.RawUI.ForegroundColor = $OriginalColor
}

function Write-Info {
    param([string]$Message)
    Write-ColoredOutput "Blue" "INFO" $Message
}

function Write-Success {
    param([string]$Message)
    Write-ColoredOutput "Green" "SUCCESS" $Message
}

function Write-Warning {
    param([string]$Message)
    Write-ColoredOutput "Yellow" "WARNING" $Message
}

function Write-Error {
    param([string]$Message)
    Write-ColoredOutput "Red" "ERROR" $Message
}

function Show-Help {
    Write-Host "SPH Particle Simulator Build Script (Windows)"
    Write-Host ""
    Write-Host "Usage: .\build.ps1 [OPTIONS]"
    Write-Host ""
    Write-Host "Options:"
    Write-Host "  -Debug          Build in Debug mode"
    Write-Host "  -Release        Build in Release mode (default)"
    Write-Host "  -NoOpenMP       Disable OpenMP parallelization"
    Write-Host "  -CUDA           Enable CUDA GPU acceleration"
    Write-Host "  -NoPython       Disable Python bindings"
    Write-Host "  -NoExamples     Don't build examples"
    Write-Host "  -NoBenchmarks   Don't build benchmarks"
    Write-Host "  -Tests          Build unit tests"
    Write-Host "  -Verbose        Verbose output"
    Write-Host "  -Clean          Clean build directory first"
    Write-Host "  -Help           Show this help message"
    Write-Host ""
    Write-Host "Examples:"
    Write-Host "  .\build.ps1                    # Standard release build"
    Write-Host "  .\build.ps1 -Debug -Tests      # Debug build with tests"
    Write-Host "  .\build.ps1 -CUDA              # Release build with CUDA"
    Write-Host "  .\build.ps1 -Clean             # Clean and rebuild"
    exit 0
}

function Test-Dependencies {
    Write-Info "Checking system dependencies..."

    $missingTools = @()

    # Check for CMake
    if (!(Get-Command cmake -ErrorAction SilentlyContinue)) {
        $missingTools += "cmake"
    }

    # Check for Visual Studio or MSVC
    $vsPath = "${env:ProgramFiles(x86)}\Microsoft Visual Studio"
    if (!(Test-Path $vsPath)) {
        $missingTools += "Visual Studio (MSVC compiler)"
    }

    # Check for vcpkg or dependencies
    if (!(Get-Command vcpkg -ErrorAction SilentlyContinue)) {
        Write-Warning "vcpkg not found - you'll need to install GLFW, GLM manually"
    }

    if ($missingTools.Count -gt 0) {
        Write-Error "Missing dependencies: $($missingTools -join ', ')"
        Write-Info "Please install them and run this script again."
        Write-Info "Recommended: Install Visual Studio with C++ workload and use vcpkg for dependencies"
        exit 1
    }

    Write-Success "All basic dependencies found!"
}

function Clear-BuildDirectory {
    if ($CleanBuild -or !(Test-Path $BuildDir)) {
        Write-Info "Cleaning build directory..."
        if (Test-Path $BuildDir) {
            Remove-Item -Recurse -Force $BuildDir
        }
        New-Item -ItemType Directory -Path $BuildDir | Out-Null
        Write-Success "Build directory cleaned"
    }
}

function Invoke-CMakeConfigure {
    Write-Info "Configuring CMake..."

    $cmakeArgs = @(
        "-DCMAKE_BUILD_TYPE=$BuildType",
        "-DUSE_OPENMP=$UseOpenMP",
        "-DUSE_CUDA=$UseCUDA",
        "-DBUILD_PYTHON_BINDINGS=$BuildPython",
        "-DBUILD_EXAMPLES=$BuildExamples",
        "-DBUILD_BENCHMARKS=$BuildBenchmarks",
        "-DBUILD_TESTS=$BuildTests"
    )

    if ($VerboseMode -eq "ON") {
        $cmakeArgs += "-DCMAKE_VERBOSE_MAKEFILE=ON"
    }

    Push-Location $BuildDir

    try {
        & cmake .. @cmakeArgs
        if ($LASTEXITCODE -eq 0) {
            Write-Success "CMake configuration completed"
        } else {
            Write-Error "CMake configuration failed"
            exit 1
        }
    } finally {
        Pop-Location
    }
}

function Invoke-Build {
    Write-Info "Building project..."

    Push-Location $BuildDir

    try {
        $buildArgs = @("--build", ".", "--config", $BuildType, "--parallel")
        if ($VerboseMode -eq "ON") {
            $buildArgs += "--verbose"
        }

        & cmake @buildArgs

        if ($LASTEXITCODE -eq 0) {
            Write-Success "Build completed successfully!"
        } else {
            Write-Error "Build failed"
            exit 1
        }
    } finally {
        Pop-Location
    }
}

function Invoke-Tests {
    if ($BuildTests -eq "ON") {
        Write-Info "Running tests..."

        Push-Location $BuildDir

        try {
            & ctest --output-on-failure --build-config $BuildType
            if ($LASTEXITCODE -eq 0) {
                Write-Success "All tests passed!"
            } else {
                Write-Warning "Some tests failed"
            }
        } catch {
            Write-Warning "ctest not available or failed to run"
        } finally {
            Pop-Location
        }
    }
}

function Invoke-Benchmarks {
    if ($BuildBenchmarks -eq "ON") {
        Write-Info "Running benchmarks..."

        $benchmarkPath = Join-Path $BuildDir "benchmarks\performance_test.exe"
        if (Test-Path $benchmarkPath) {
            & $benchmarkPath 1000 5000 10000 25000
        } else {
            Write-Warning "Benchmark executable not found"
        }
    }
}

function Show-Summary {
    Write-Success "Build Summary:"
    Write-Host "  Build Type: $BuildType"
    Write-Host "  OpenMP: $UseOpenMP"
    Write-Host "  CUDA: $UseCUDA"
    Write-Host "  Python Bindings: $BuildPython"
    Write-Host "  Examples: $BuildExamples"
    Write-Host "  Benchmarks: $BuildBenchmarks"
    Write-Host "  Tests: $BuildTests"
    Write-Host ""

    Write-Info "Executables created in: $BuildDir\"
    $exePath = Join-Path $BuildDir "Release\sph_simulator.exe"
    if (Test-Path $exePath) {
        Write-Host "  - sph_simulator.exe (main executable)"
    }

    $examplesPath = Join-Path $BuildDir "examples"
    if (Test-Path $examplesPath) {
        Write-Host "  - examples\ (simulation examples)"
    }

    $benchmarksPath = Join-Path $BuildDir "benchmarks"
    if (Test-Path $benchmarksPath) {
        Write-Host "  - benchmarks\ (performance tests)"
    }

    if ($BuildPython -eq "ON") {
        Write-Host "  - python\ (Python bindings)"
    }

    Write-Host ""
    Write-Info "To run the simulator:"
    Write-Host "  cd $BuildDir\Release"
    Write-Host "  .\sph_simulator.exe"
    Write-Host ""
    Write-Info "To run examples:"
    Write-Host "  .\examples\dam_break.exe 10000"
}

# Main script logic
if ($Help) {
    Show-Help
}

Write-Info "SPH Particle Simulator Build Script (Windows)"
Write-Info "PowerShell Version: $($PSVersionTable.PSVersion)"

# Parse arguments
if ($Debug) { $BuildType = "Debug" }
if ($NoOpenMP) { $UseOpenMP = "OFF" }
if ($CUDA) { $UseCUDA = "ON" }
if ($NoPython) { $BuildPython = "OFF" }
if ($NoExamples) { $BuildExamples = "OFF" }
if ($NoBenchmarks) { $BuildBenchmarks = "OFF" }
if ($Tests) { $BuildTests = "ON" }
if ($Verbose) { $VerboseMode = "ON" }
if ($Clean) { $CleanBuild = $true }

Test-Dependencies
Clear-BuildDirectory
Invoke-CMakeConfigure
Invoke-Build
Invoke-Tests
Invoke-Benchmarks
Show-Summary

Write-Success "Build process completed successfully! 🚀"

