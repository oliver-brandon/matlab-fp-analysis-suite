data.epocs.cNoRewA.name = 'cNoRewA'
data.epocs.cNoRewA.onset = onset
data.epocs.cNoRewA.offset = onset + 1
data.epocs.cNoRewA.data = ones(height(onset))

data.epocs.iNoRewA.name = 'iNoRewA'
data.epocs.iNoRewA.onset = onset
data.epocs.iNoRewA.offset = onset + 1
data.epocs.iNoRewA.data = ones(height(onset))

data.epocs.cRewA.name = 'cRewA'
data.epocs.cRewA.onset = onset
data.epocs.cRewA.offset = onset + 1
data.epocs.cRewA.data = ones(height(onset))

data.epocs.iRewA.name = 'iRewA'
data.epocs.iRewA.onset = onset
data.epocs.iRewA.offset = onset + 1
data.epocs.iRewA.data = ones(height(onset))

save("E:\Google Drive\optomouse-prime\mats-unprocessed\1035F_Rev1Sham.mat","data")