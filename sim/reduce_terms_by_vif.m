function [keepMask, removalLog] = reduce_terms_by_vif(A, exponents, threshold)
%REDUCE_TERMS_BY_VIF すべての項のVIFがしきい値未満になるまで、高次項を1つずつ削除する。
%   削除する項の選び方
%     1) VIFがしきい値以上の高次項があれば、その中で次数が最も高い項（同じ次数ならVIF最大）
%     2) しきい値以上が1次項だけなら、高次項の中でVIF最大の項（1次項と相関している項）
%   VIF最大の項を単純に消すと y^3 を消して y^6 を残すような結果になるため、
%   低次の項を優先して残す（多項式モデルの階層性を保つ）ようにしている。
%   切片と1次項（並進・倍率・回転）は物理的に必須なので削除しない。
%   keepMask   : 残す項なら true
%   removalLog : [削除した項の番号, 削除直前のVIF]（削除した順）
keepMask = true(size(exponents, 1), 1);
isProtected = sum(exponents, 2) <= 1;
removalLog = zeros(0, 2);
while true
    keptIndex = find(keepMask);
    vif = compute_vif(A(:, keptIndex), exponents(keptIndex, :));
    if max(vif, [], 'omitnan') < threshold
        break;
    end
    removableVif = vif;
    removableVif(isProtected(keptIndex)) = -Inf;
    isOverThreshold = removableVif >= threshold;
    if any(isOverThreshold)
        candidate = find(isOverThreshold);
        candidateOrder = sum(exponents(keptIndex(candidate), :), 2);
        candidate = candidate(candidateOrder == max(candidateOrder));
        [maxVif, bestInCandidate] = max(removableVif(candidate));
        position = candidate(bestInCandidate);
    else
        [maxVif, position] = max(removableVif);
    end
    if ~isfinite(maxVif)
        break;  % 削除できる高次項が残っていない
    end
    removalLog(end + 1, :) = [keptIndex(position), maxVif]; %#ok<AGROW>
    keepMask(keptIndex(position)) = false;
end
end
