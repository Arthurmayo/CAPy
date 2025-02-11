import os

def convert_awd_to_csv(awd_file_path, csv_file_path):
    """
    Convert a .awd file to a cleaned .csv file for use with CAPy.
    
    The conversion does the following:
      - Keeps the first 4 lines as is.
      - Clears (empties) lines 5, 6, and 7.
      - For every subsequent line:
          - Splits the line by whitespace and takes only the first column.
          - If that value is negative, the entire row is skipped (deleted).
    
    Parameters:
      awd_file_path (str): Path to the input .awd file.
      csv_file_path (str): Path to the output .csv file.
    
    Returns:
      csv_file_path (str): The output CSV file path if conversion is successful.
    
    Raises:
      IOError: If reading or writing the file fails.
    """
    try:
        with open(awd_file_path, 'r') as infile:
            lines = infile.readlines()
    except Exception as e:
        raise IOError(f"Error reading AWD file: {e}")
    
    cleaned_lines = []
    for i, line in enumerate(lines):
        if i < 4:
            # Keep the first 4 lines as is.
            cleaned_lines.append(line.rstrip() + '\n')
        elif 4 <= i < 7:
            # Clear lines 5, 6, and 7 by writing an empty line.
            cleaned_lines.append('\n')
        else:
            # For every subsequent line, split by whitespace and take the first column.
            parts = line.strip().split()
            if parts:
                try:
                    value = float(parts[0])
                except ValueError:
                    # If conversion fails, skip the row.
                    continue
                # Skip the row if the value is negative.
                if value < 0:
                    continue
                # Otherwise, write only the first value.
                cleaned_lines.append(f"{parts[0]}\n")
            else:
                cleaned_lines.append('\n')
    
    try:
        with open(csv_file_path, 'w') as outfile:
            outfile.writelines(cleaned_lines)
    except Exception as e:
        raise IOError(f"Error writing CSV file: {e}")
    
    return csv_file_path

if __name__ == "__main__":
    # Example usage from the command line:
    # python awd_converter.py input.awd output.csv
    import sys
    if len(sys.argv) != 3:
        print("Usage: python awd_converter.py input.awd output.csv")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        try:
            convert_awd_to_csv(input_file, output_file)
            print(f"Conversion successful: {output_file}")
        except Exception as e:
            print(f"Conversion failed: {e}")