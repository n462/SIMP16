# Find first value that is below threshold
def thres(values, cutoff):
    trunc = len(values)
    for value in values:
        if value < cutoff:
            trunc = values.index(value)
            break
    return trunc

# Return the trimming position 
def trunc(F, R, paired):
    import pandas as pd
    dataF = pd.read_table(F)
    qmin = 25# quality threshold
    rmin = 0.75 # read count threshold
    
    # select row with score values
    scoreF = dataF.iloc[0][1:].values.tolist()
    # select row with count numbers
    countF = dataF.iloc[1][1:].values.tolist()
    # extract maximum read count
    max_readsF = int(max(countF) * rmin)
    # find quality threshold
    trunc_scoreF = thres(scoreF, qmin)
    # find read threshold
    trunc_countF = thres(countF, max_readsF)
    # final threshold
    truncF = min(trunc_scoreF, trunc_countF)
    
    if paired:
        dataR = pd.read_table(R)
        scoreR = dataR.iloc[0][1:].values.tolist()
        countR = dataR.iloc[1][1:].values.tolist()
        max_readsR = int(max(countR) * rmin)
        trunc_scoreR = thres(scoreR, qmin)
        trunc_countR = thres(countR, max_readsR)
        truncR = min(trunc_scoreR, trunc_countR)
        return(truncF,truncR)
    else:
        return(truncF,)
