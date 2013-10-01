library(zinba)

generateAlignability(
   outdir = '~/data/esc/faire/out/',
   mapdir = '~/data/littman/FAIRE_seq/map50_mm9/',
   athresh = 1,
   extension = 36,
   twoBitFile = '~/data/littman/FAIRE_seq/mm9.2bit'
   )

basealigncount(
  inputfile='Sample_lane6.E14.sorted.bed', #mapped sample reads
  outputfile= '~/data/littman/FAIRE_seq/out/', # output path
  extension= 36, #average fragment library length
  filetype='bed', #either "bed", "bowtie", or "tagAlign"
  twoBitFile='~/data/littman/FAIRE_seq/mm9.2bit', #path to downloaded genome build file in .2bit format
)

run.zinba(
  seq='Sample_lane6.E14.sorted.bed',
  input='none',
  filetype='bed',
  twoBit='mm9.2bit',
  winSize=250,
  offset=125,
  extension=200,
  basecountfile='Sample_lane6.E14.sorted.basecount',
  align='out/',
  selectmodel=T,
  selectchr='chr20',
  selecttype='dirty',
  selectcovs=c('gcPerc','align_perc','exp_cnvwin_log'),
  interaction=T,
  threshold=0.05,
  refinepeaks=1,
  numProc=2,
  winGap=0,
  FDR=TRUE,
  outfile='zinba_out/test',
  printFullOut=1,
  method='mixture'
  )

run.zinba(
  + seq = 'SL14401.bed',
  + input = 'none',
  + filetype='bed',
  + twoBit='mm9.2bit',
  + winSize=500,
  + offset=250,
  + extension=200,
  + basecountfile='SL14401.basecount',
  + align='out/',
  + selectmodel=F,
  + formula = exp_count ~ align_perc + exp_cnvwin_log + gcPerc + align_perc:exp_cnvwin_log + align_perc:gcPerc + exp_cnvwin_log:gcPerc + align_perc:exp_cnvwin_log:gcPerc + 1,
  + formulaE = exp_count ~ align_perc + exp_cnvwin_log + gcPerc + align_perc:exp_cnvwin_log + 1,
  + formulaZ = exp_count ~ align_perc + exp_cnvwin_log + gcPerc + align_perc:exp_cnvwin_log + 1,
  + threshold=0.05,
  + refinepeaks=1,
  + numProc=2,
  + winGap=0,
  + FDR=TRUE,
  + outfile='zinba_out/rerun',
  + printFullOut=1,
  + method='mixture'
  + )
