function output = count_unique_entries(cell_array)

% concatenate all the cell arrays into one big cell array
all_entries = vertcat(cell_array{:});

% count number of unique entries
unique_entries = unique(all_entries);
num_unique_entries = length(unique_entries);

% count occurrences of each unique entry
counts = zeros(num_unique_entries, 1);
for i = 1:num_unique_entries
    counts(i) = sum(strcmp(all_entries, unique_entries{i}));
end

% sort in descending order by count
[sorted_counts, sorted_indices] = sort(counts, 'descend');
sorted_entries = unique_entries(sorted_indices);

% output results
output = table(sorted_entries, sorted_counts, 'VariableNames', {'Entry', 'Count'});

end