module Frost

include("TanTRG.jl")
include("XTRG.jl")

export TanTRG, XTRG, tantrg, xtrg

const tantrg = TanTRG.main
const xtrg = XTRG.xtrg_main
const DEFAULT_LATTICE_INFO_PATH = TanTRG.DEFAULT_LATTICE_INFO_PATH

end