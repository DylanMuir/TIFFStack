%%
vfSignal = mfRegionTraces(4, :);

vfCTrans = CalciumTransientTemplate(3.91);

nNumCorrs = numel(vfSignal)-numel(vfCTrans);

mfSignal = zeros(nNumCorrs, numel(vfCTrans));
mfCorrs = zeros(nNumCorrs, numel(vfCTrans));

for i = 1:nNumCorrs
   mfSignal(i, :) = vfSignal(i:i+numel(vfCTrans)-1);
   vfCorr = fftshift(xcov(mfSignal(i, :), vfCTrans));
   mfCorrs(i, :) = vfCorr(1:numel(vfCTrans));
end

%% Perform correlation

figure;
plot(mfCorrs(:, 1));
