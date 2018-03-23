module AlignUtils

Pkg.installed("DISSEQT")==nothing && warn("Module DISSEQT not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")
Pkg.installed("SynapseClient")==nothing && warn("Module SynapseClient not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")
Pkg.installed("SynapseTools")==nothing && warn("Module SynapseTools not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")
Pkg.installed("BamReader")==nothing && warn("Module BamReader not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")


using DataFrames
using Gadfly
using Colors
using JLD
Pkg.installed("MAT")!=nothing && using MAT
using Levenshtein
using Base.Collections
using SynapseClient
using SynapseTools
using BamReader
using DISSEQT


export
	Sample,
	find_samples,
	find_aligned,
	assign_reference!,
	align_sample!,
	align_samples!,
	askforconfirmation,
	makecleanfolder,
	getreferenceinfo,
	reference_sanity_check,
	computecodonfrequencies,
	computenucleotidefrequencies,
	uploadaligned,
	uploadswarms,
	coverageplots

include("log.jl")
include("sample.jl")
include("synapseutils.jl")
include("align.jl")
include("misc.jl")
include("consensus.jl")
include("codonfrequency.jl")
include("nucleotidefrequency.jl")
include("readcoverage.jl")

end
