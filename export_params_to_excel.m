function export_params_to_excel(params, filename)

    % Get field names
    fields = fieldnames(params);

    % Determine number of beats (length of any vector field)
    N = [];
    for f = 1:numel(fields)
        val = params.(fields{f});
        if isnumeric(val) && isvector(val) && numel(val) > 1
            N = numel(val);
            break
        end
    end

    if isempty(N)
        error('No vector fields found in params. Nothing to export per beat.');
    end

    % Build table row-by-row
    T = table();
    for f = 1:numel(fields)
        name = fields{f};
        val  = params.(name);

        if isscalar(val)
            % replicate scalar across all beats
            T.(name) = repmat(val, N, 1);
        elseif isvector(val) && numel(val) == N
            % use vector as-is
            T.(name) = val(:);
        else
            warning('Skipping field "%s" due to incompatible size.', name);
        end
    end

    % Write to Excel
    writetable(T, filename);
    fprintf('Exported %d rows to %s\n', N, filename);

end