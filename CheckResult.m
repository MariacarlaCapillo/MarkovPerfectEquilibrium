function check = CheckResult(message)
% Check if solver worked by analyzing its message

containsMessage = contains(message, 'Equation solved');
containsMessage_stalled = contains(message, 'Equation solved, fsolve stalled');
% Set check=1 if the solver didn't work properly, i.e. if the following
% messages are displayed.
if containsMessage == 0 || containsMessage_stalled == 1
    check = 1;
else
    check = 0;
end

end

