from scipy.stats import ttest_ind
import pandas as pd

# Load data
df = pd.read_csv(r'D:\kiel\theses\master_theses_lgodbersen\data\PuG\participants.tsv', sep='\t')

# Convert 'age' to numeric if necessary
df['age'] = pd.to_numeric(df['age'], errors='coerce')

# exclude the following IDs
ids_to_exclude = ['sub-AN20SC', 'sub-BI19KI', 'sub-CA26NM', 'sub-EL16KI', 'sub-HE25KI', 'sub-IM31KI', 'sub-JA04NE', 'sub-KA05KI', 'sub-KA08LU', 'sub-MA26WO', 'sub-MO09MU', 'sub-UT16FL', 'sub-WI17PE','sub-AN14NE', 'sub-HI22NE', 'sub-RE03IT', 'sub-DO06CL']
df = df[~df['participant_id'].isin(ids_to_exclude)]

# Split the DataFrame
df_withPCS = df[df['group'] == 'withPCS']
df_withoutPCS = df[df['group'] == 'withoutPCS']

# Function to find the closest match based on age
def find_closest_match(target_age, df_candidates):
    age_diff = (df_candidates['age'] - target_age).abs()
    closest_index = age_diff.idxmin()
    return closest_index

matched_withPCS_indices = []
for index, row in df_withoutPCS.iterrows():
    # Find the closest match
    closest_index = find_closest_match(row['age'], df_withPCS.drop(matched_withPCS_indices))
    
    # Append the index of the match to the list
    matched_withPCS_indices.append(closest_index)

# Extract matched 'withPCS' participants
matched_withPCS = df_withPCS.loc[matched_withPCS_indices]

# Perform the t-test on ages
t_statistic, p_value = ttest_ind(matched_withPCS['age'], df_withoutPCS['age'], nan_policy='omit')

# Print results
print("T-value: ", t_statistic)
print("P-value: ", p_value)
print("Group size ('withPCS'): ", matched_withPCS['age'].count())
print("Group size ('withoutPCS'): ", df_withoutPCS['age'].count())

# Save to TSV
combined_df = pd.concat([matched_withPCS, df_withoutPCS])
combined_df.to_csv(r'D:\kiel\theses\master_theses_lgodbersen\data\PuG\matched_participants_conn.tsv', sep='\t', index=False)
