# Skye Goetz (CalPoly) 10/22/2024

using FilePathsBase
using LaTeXStrings
using DataFrames
using PyCall
using Plots
using YAML
using CSV

function SmilesTo3DCoordinates(SmilesFormula::String)::Vector{String}
    # Converts SmilesFormula to 3D coordinates with a PyCall to RDKit

    Chem::Any = pyimport("rdkit.Chem")
    AllChem::Any = pyimport("rdkit.Chem.AllChem")    
    
    # Builds RDKit structure and conformer
    BasicSmilesStructure::Any = ""
    try
        BasicSmilesStructure = Chem.MolFromSmiles(SmilesFormula)
    catch e
        SmilesErrorMessage::String = "Error: Cannot Recognize Smiles -> $(e)"
        println(SmilesErrorMessage)
        exit(1)
    end
    SmilesStructure::Any = ""
    try
        SmilesStructure = Chem.AddHs(BasicSmilesStructure, explicitOnly=true)
    catch e
        HydrogensErrorMessage::String = "Error: Cannot Add Hydrogens -> $(e)"
        println(HydrogensErrorMessage)
        exit(1)
    end
    AllChem.EmbedMolecule(SmilesStructure)
    AllChem.UFFOptimizeMolecule(SmilesStructure)
    SmilesConformer::Any = SmilesStructure.GetConformer()
    
    # Extracts coordinates from conformer
    XYZCoordinates::Vector{String} = []
    for Atom::Any in SmilesStructure.GetAtoms()
        Position::Any = SmilesConformer.GetAtomPosition(Atom.GetIdx())
        AtomicSymbol::String = Atom.GetSymbol()
        X::Float64, Y::Float64, Z::Float64 = Position.x, Position.y, Position.z
        AtomXYZ::String = "$(AtomicSymbol)\t$(X)\t$(Y)\t$(Z)\n"
        push!(XYZCoordinates, AtomXYZ)
    end

    return XYZCoordinates
end

function QueryWriter(Type::String, Coordinates::Vector{String}, MoleculeSubConfig::Dict{Any, Any}, CommandSubConfig::Dict{Any, Any}, IterationNumber::Int64)::String
    # Creates and writes coordinates and commands to an Orca input file

    Functional = CommandSubConfig["Functional"]
    CleanedSmiles::String = replace(MoleculeSubConfig["SmilesFormula"], r"[^a-zA-Z]" => "")
    PathToOrcaQuery::String = "Orca.$(IterationNumber).$(Type).$(CleanedSmiles).$(Functional).inp"
    BasisSet::String = CommandSubConfig["BasisSet"]
    SolventModel::String = CommandSubConfig["SolventModel"]
    Solvent::String = CommandSubConfig["SolventToUse"]

    # Creates case-specific type blocks
    Operation::String = ""
    TypeBlock::String = ""
    
    if Type == "TightOpt"
        Operation = "TightOpt"
        TypeBlock = "\n"
    elseif Type == "ScanAngle"
        Operation = "TightOpt"
        ScanParameters::Vector{String} = []
        ScanType::String = CommandSubConfig["ScanType"]
        push!(ScanParameters, ScanType)
        for Angle in CommandSubConfig["Angles"]
            push!(ScanParameters, string(Angle))
        end
        push!(ScanParameters, "=")
        for RangeLimit in CommandSubConfig["Advanced"]["ScanRange"]
            push!(ScanParameters, "$RangeLimit,")
        end     
        NumberOfScans::Any = CommandSubConfig["Advanced"]["NumberOfScans"]
        push!(ScanParameters, string(NumberOfScans))
        ScanMode::String = join(ScanParameters, " ")
        TypeBlock = "\n%geom Scan\n$(ScanMode)\nend\nend\n\n"
    elseif Type == "TDDFT"
        Operation = ""
        NumRoots::Int64 = CommandSubConfig["Advanced"]["Nroots"]
        Maxdim::Int64 = CommandSubConfig["Advanced"]["Maxdim"]
        TypeBlock = "\n%tddft\nnroots $(NumRoots)\nmaxdim $(Maxdim)\nend\n\n"
    end
    
    # Inserts query parameters to an array
    Charge::Int64 = MoleculeSubConfig["MolecularCharge"]
    Multiplicity::Int64 = MoleculeSubConfig["MolecularMultiplicity"]
    PreCoordinates::String = "* xyz $(string(Charge)) $(string(Multiplicity))\n"
    insert!(Coordinates, 1, PreCoordinates)
    insert!(Coordinates, 1, TypeBlock)
    CommandLine::String = "! $(Functional) $(BasisSet) $(Operation) $(SolventModel)($(Solvent))\n"
    insert!(Coordinates, 1, CommandLine)
    push!(Coordinates, "*\n")
    
    # Writes Orca 6 query from the array
    open(PathToOrcaQuery, "w") do Query::IO
        for Line in Coordinates
            write(Query, Line)
        end
    end

    return PathToOrcaQuery
end

function ReplaceExtension(Path::String, NewExtension::String)::String
    # Replaces Extensions

    BasePath::String, _ = splitext(Path)
    return BasePath * "." * NewExtension
end

function RunQuery(Mode::String, Config::Dict{String, Any}, Coordinates::Vector{String}, Molecule::Any, IterationNumber)::String
    # Executes Orca Query and Captures Output

    OutExtension::String = "out"
    LogExtension::String = "log"

    PathToQuery::String = QueryWriter(Mode, Coordinates, Molecule, Config[Mode], IterationNumber)
    PathToOutput::String = ReplaceExtension(PathToQuery, OutExtension)
    
    # Runs Orca if Output Doesn't Already Exist
    if !isfile(PathToOutput)
    PathToLog::String = ReplaceExtension(PathToQuery, LogExtension)
    OrcaCommand::Cmd = `orca $(PathToQuery)`
        open(PathToOutput, "w") do Output::IO
            open(PathToLog, "w") do ErrorLog::IO
                run(pipeline(OrcaCommand, stdout=Output, stderr=ErrorLog))
                Stderr::String = read(ErrorLog, String)
                if !isempty(Stderr)
                    StderrMessage::String = "$(Mode) Error: $(Stderr)"
                    println(StderrMessage)
                    exit(1)
                end
            end
        end
    end

    return PathToOutput
end

function CoordinateRegex(S::String)::String
    # Regex For Extracting Coordinates From Orca

    return replace(S, r"^\s+|\s+$" => "")
end

function TightOptParser(PathToTightOptOutput::String)::Vector{String}
    # Parses TightOpt Output From Orca

    TightOptCoordinates::Vector = []
    Count::Int64 = -1

    XYZExtension::String = "xyz"
    XYZFile::String = ReplaceExtension(PathToTightOptOutput, XYZExtension)
    
    open(XYZFile) do Output::IO
        for Line::String in eachline(Output)
            Count += 1
            
            # Skips First Two Lines In The File
            if Count > 1
                push!(TightOptCoordinates, CoordinateRegex(Line) * "\n")
            end
        end
    end

    return TightOptCoordinates
end

function ScanAngleParser(PathToScanAngleOutput::String, Coordinates::Vector{String})::Dict{Int64, Vector{String}}
    # Parses ScanAngle Output From Orca

    # Initialize the dictionary
    ScanAngleCoordinates::Dict{Int64, Vector{String}} = Dict()

    StringToStartBlock::String = "FINAL ENERGY EVALUATION AT THE STATIONARY POINT"
    StringToStartCoordinates::String = "CARTESIAN COORDINATES (ANGSTROEM)"

    IsNewBlock::Bool = false
    BlockCounter::Int64 = 0
    IsCoordinates::Bool = false
    CoordinateCounter::Int64 = 0
    CoordinateList::Vector{String} = []

    open(PathToScanAngleOutput) do Output::IO
        for Line::String in eachline(Output)
            
            # Checks For Convergence at A Scan
            if occursin(StringToStartBlock, Line)
                IsNewBlock = true
                IsCoordinates = false
                BlockCounter += 1
            end

            # Parses Coordinates ONLY at each successful scan
            if IsNewBlock
                if occursin(StringToStartCoordinates, Line)
                    IsCoordinates = true
                    # To Skip Lines At The Beginning Of Each Coordinate Block
                    CoordinateCounter = -2
                    CoordinateList = []
                end

                if IsCoordinates
                    CoordinateCounter += 1
                    if length(Coordinates) >= CoordinateCounter && CoordinateCounter > 0
                        push!(CoordinateList, CoordinateRegex(Line) * "\n")
                    elseif CoordinateCounter > length(Coordinates)
                        IsNewBlock = false
                        IsCoordinates = false
                        
                        # Add CoordinateList to the dictionary with BlockCounter as the key
                        if haskey(ScanAngleCoordinates, BlockCounter)
                            push!(ScanAngleCoordinates[BlockCounter], CoordinateList...)
                        else
                            ScanAngleCoordinates[BlockCounter] = CoordinateList
                        end
                    end
                end
            end
        end
    end

    return ScanAngleCoordinates
end

function TDDFTPaser(PathToTDDFTOutput::String)::Tuple{Vector{Float64}, Vector{Float64}}
    # Parses TDDFT Output From Orca

    StringToBeginSpectra::String = "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS"
    StringToEndSpectra::String = "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"

    Fosc::Vector{Float64} = []
    Wavelength::Vector{Float64} = []
    IsBlock::Bool = false
    open(PathToTDDFTOutput) do Output::IO
        for Line::String in eachline(Output)
            # Checks For Coordinate Block
            if occursin(StringToBeginSpectra, Line)
                IsBlock = true
            elseif occursin(StringToEndSpectra, Line)
                IsBlock = false
            end
            if IsBlock
                TableRow::Vector{Any} = split(Line)
                # We have to tryparse bc everything's read as a String
                if length(TableRow) >= 6 && tryparse(Float64, TableRow[6]) !== nothing 
                    push!(Fosc, parse(Float64, TableRow[7]))
                    push!(Wavelength, parse(Float64, TableRow[6]))                    
                end
            end
        end
    end

    return Fosc, Wavelength
end

function TTDFTToCSV(Wavelength::Vector{Float64}, Fosc::Vector{Float64}, IterationNumber::Int64)::Nothing
    # Saves TDDFT Values Used In Spectra Graphs 

    PathToTDDFTCSVs::String = "TDDFTGraphData"
    mkpath(PathToTDDFTCSVs)

    CSVPath::String = "TDDFTGraphData.$(IterationNumber).csv"
    Frame::DataFrames.DataFrame = DataFrame(Wavelength = Wavelength, Fosc = Fosc)
    CSV.write(CSVPath, Frame)
    mv(CSVPath, joinpath(PathToTDDFTCSVs, CSVPath), force=true)

    return nothing
end

function GraphTDDFT(X::Vector{Float64}, Y::Vector{Float64}, IterationNumber::Int64)::Nothing
    # Graphs Wavelength and Fosc From TDDFT Calculations, Saves Graph Values

    PathToTDDFTGraphs::String = "TDDFTGraphs"
    mkpath(PathToTDDFTGraphs)

    bar(X, Y, width=0.25, color=:blue, legend=false, stroke=:blue)
    plot!(title="Absorbance Spectra $(IterationNumber)", xlabel=L"Wavelength (nm)", ylabel=L"Fosc\ (au^{2})", grid=:y, gridalpha=0.7, gridlinestyle=:dash)    

    PathToFigure::String = "TDDFTSpectra.$(IterationNumber).png"
    savefig(PathToFigure)
    mv(PathToFigure, joinpath(PathToTDDFTGraphs, PathToFigure), force=true)

    TTDFTToCSV(X, Y, IterationNumber)

    return nothing
end

function XYZWriter(Mode::String, Coordinates::Vector{String}, IterationNumber::Int64)::Nothing
    # Writes Coordinates to an XYZ File

    XYZDirectory::String = "XZY"
    mkpath(XYZDirectory)

    XYZFilePath::String = "$(Mode).$(IterationNumber).xyz"
    open(XYZFilePath, "w") do XYZ::IO
        for Coordinate::String in Coordinates
            write(XYZ, Coordinate)
        end
    end
    mv(XYZFilePath, joinpath(XYZDirectory, XYZFilePath), force=true)

    return nothing
end

function ProcessMolecule(Molecule::Any, Config::Dict{String, Any})::Nothing
    # Block to write and process Orca commands 

    BaseIterationNumber::Int64 = 1
    XYZCoordinates::Vector{String} = SmilesTo3DCoordinates(Molecule["SmilesFormula"])
    XYZWriter("Initial", XYZCoordinates, BaseIterationNumber)

    # TightOpt
    PathToTightOptOutput::String = RunQuery("TightOpt", Config, copy(XYZCoordinates), Molecule, BaseIterationNumber)
    TightOptCoordinates::Vector{String} = TightOptParser(PathToTightOptOutput)
    XYZWriter("TightOpt", copy(TightOptCoordinates), BaseIterationNumber)

    PathToTDDFTOutput::String = RunQuery("TDDFT", Config, copy(TightOptCoordinates), Molecule, BaseIterationNumber)
    InitialFosc::Vector{Float64}, InitialWavelength::Vector{Float64} = TDDFTPaser(PathToTDDFTOutput)
    GraphTDDFT(InitialWavelength, InitialFosc, BaseIterationNumber)

    # ScanAngle
    PathToScanAngleOutput::String = RunQuery("ScanAngle", Config, copy(TightOptCoordinates), Molecule, BaseIterationNumber)
    ScanAngleCoordinates::Dict{Int64, Vector{String}} = ScanAngleParser(PathToScanAngleOutput, copy(TightOptCoordinates))

    # TDDFT
    TDDFTCounter::Int64 = 0
    for (Step::Int64, ScanCoordinates::Vector{String}) in ScanAngleCoordinates
        TDDFTCounter += 1
        XYZWriter("ScanAngle", copy(ScanCoordinates), Step)
        if TDDFTCounter == Config["AdvancedSettings"]["NumberOfGeometryScansPerTDDFTCalculation"]
            PathToTDDFTOutput = RunQuery("TDDFT", Config, copy(ScanCoordinates), Molecule, Step)
            Fosc::Vector{Float64}, Wavelength::Vector{Float64} = TDDFTPaser(PathToTDDFTOutput)
            GraphTDDFT(Wavelength, Fosc, Step)
            TDDFTCounter = 0
        end
    end
    
    return nothing
end

function MainBlock(PathToYAML::String)::Nothing
    # Main block for OrcaFunctionalHub

    Config::Dict{String, Any} = Dict()
    try
        Config = YAML.load_file(PathToYAML)
    catch e
        YAMLErrorMessage::String = "Error: Cannot Load YAML -> $(e)"
        println(YAMLErrorMessage)
        exit(1)
    end

    for Molecule::Any in Config["Molecules"]
        ProcessMolecule(Molecule, Config)
    end

    return nothing
end

# Execution block
if length(ARGS) != 1
    UsageMessage::String = "Usage: julia OrcaFunctionalHub.jl <PathToYAML>"
    println(UsageMessage)
    exit(1)
else
    PathToYAML::String = ARGS[1]
    MainBlock(PathToYAML)
end