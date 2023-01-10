% trials is channels x time x trials
function trials = offsetTrials(trials,zeroT)

% Offset
for i1 = 1:size(trials,1)
    for i2 = 1:size(trials,3)
        offset = trials(i1,zeroT,i2);
        trials(i1,:,i2) = trials(i1,:,i2) - offset;
    end
end

end