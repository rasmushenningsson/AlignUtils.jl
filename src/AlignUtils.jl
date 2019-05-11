module AlignUtils

using Pkg

haskey(Pkg.installed(),"DISSEQT") || warn("Module DISSEQT not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")
haskey(Pkg.installed(),"SynapseClient") || warn("Module SynapseClient not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")
haskey(Pkg.installed(),"SynapseTools") || warn("Module SynapseTools not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")
haskey(Pkg.installed(),"BamReader") || warn("Module BamReader not installed. Please refer to AlignUtils installation instructions at https://github.com/rasmushenningsson/AlignUtils.jl")


using Distributed
using DataFrames
using Gadfly
using Colors
using JLD
using Levenshtein
haskey(Pkg.installed(),"SynapseClient") && using SynapseClient
haskey(Pkg.installed(),"SynapseTools") && using SynapseTools
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
