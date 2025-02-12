function consistency_check_input(S, structureType)

correctFields = get_correct_fields(structureType);


userFields = fieldnames(S);

missingFields = setdiff(correctFields, userFields);

extraFields = setdiff(userFields, correctFields);

% Display the results
if ~isempty(missingFields) || ~isempty(extraFields)
    if ~isempty(missingFields)

        error(['Missing fields ' missingFields{:} ' for the ' structureType ' input structure!']);
    end
    if ~isempty(extraFields)
        error(['Missing fields ' extraFields{:} ' for the ' structureType ' input structure!']);
    end
end


end

