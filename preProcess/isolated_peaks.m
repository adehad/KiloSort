function [row, col, mu] = isolated_peaks(S1, loc_range, long_range, Th, peakLoc, postPeakSamples, maxPeakNum)
% loc_range = [3  1];
% long_range = [30  6];
smin = my_min(S1, loc_range, [1 2]);
peaks = single(S1<smin+1e-3 & S1<Th);

sum_peaks = my_sum(peaks, long_range, [1 2]);
peaks = peaks .* (sum_peaks<(maxPeakNum+0.2)).* S1;
% +0.2, because <1.2 is equivalent to <=1

peaks([1:peakLoc end-postPeakSamples:end], :) = 0;

[row, col, mu] = find(peaks);
% figure; hold on; plot(sq(gather(S1(:,1)))); plot(gather(row),sq(gather(S1(gather(row),1))),'*'); plot(sq(gather(smin(:,1)))); plot([1,size(S1,1)], [Th, Th])
    % plots trace, detected spike events, 'smoothed' spike areas and threshold