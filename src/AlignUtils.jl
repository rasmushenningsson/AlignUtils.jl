module AlignUtils

using DataFrames
using Gadfly
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
