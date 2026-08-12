#
# @Author  : zhuy
# @Date    : 2026/1/10 17:32
# @Desc    : Preset TSV file filters
#

# --- Define flexible filter factories that use column names. ---

def column_not_nan(col_name):
    """Exclude rows where the specified column is NaN."""
    return lambda df: df[col_name].notna()

def column_gt(col_name, value):
    """Filter values greater than the given threshold."""
    return lambda df: df[col_name] > value

def column_contains(col_name, keyword):
    """Filter strings containing the given keyword."""
    return lambda df: df[col_name].str.contains(keyword, na=False)

def column_not_contains(col_name, keyword):
    """Filter strings that do not contain the given keyword."""
    return lambda df: ~df[col_name].str.contains(keyword, na=False)

def column_in(col_name, options_list):
    """Filter values included in the given list."""
    return lambda df: df[col_name].isin(options_list)

def select_columns(columns):
    """Keep only the specified columns."""
    return lambda df: df[columns]

def drop_columns(columns):
    """Drop the specified columns."""
    return lambda df: df.drop(columns=columns, errors='ignore')

def regex_match(col_name, pattern):
    """Use a regular expression for advanced matching."""
    return lambda df: df[col_name].str.contains(pattern, regex=True, na=False)

def fast_assign(**kwargs):
    """
    Use Pandas assign to add multiple dynamic columns at once.
    Usage: fast_assign(Total=lambda x: x.a + x.b, Tax=lambda x: x.Total * 0.1)
    """
    return lambda df: df.assign(**kwargs)
