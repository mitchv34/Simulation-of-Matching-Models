% this function updates the value of a given parameter

function params = update_params(params, field, new_value)
    params = setfield(params, field, new_value);
end