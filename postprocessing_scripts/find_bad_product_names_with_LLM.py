#!/usr/bin/env python3

import argparse
import requests

# Note: LLAMA does a very poor job on this, GPT-4o is much better. You have to buy credit and generate a token.
def message_llm (message_user,                 
                model_source = 'GPT',
                message_system="You are a professional English linguist and a biochemist specialized on protein names. You will return an exhaustive list of product names that contain nouns in plural or repetitions of words. You know that subunit or complex are not indicators for plural. Your responses will not contain enumerations or arguments why plurals are contains. You simply return a newline-separated list of product names. Here are examples of terms that do not contain plurals: 20 kDA chaperonin chloroplastic, 20S proteasome subunit, acyltransferase. You will not return unsuspicious product names.",
                api_key = '', # insert token
               ):
    """
    This code is derived from a code snipped provided by Lars Gabriel at University
    of Greifswald.

    Sends a message to a Language Model (LLM) via an API and retrieves the response.

    This function supports communication with either the GPT3.5 Turbo or the 
    LLAMA model from University Greifswald, depending on the specified `model_source`. 

    Parameters:
    -----------
    message_user : str
        The message content from the user that will be sent to the LLM.
    
    model_source : str, optional
        The source model to use for the request, either 'GPT' or 'LLAMA'. 
        Default is 'LLAMA'.
    
    message_system : str, optional
        The system message to set the context for the conversation with the LLM.
        Default is "You are a matchmaker for scientists.".
    
    api_key : str, optional
        The API key used for authentication when making the request. 
        Default is an empty string, which requires the user to provide a valid API key.

    Raises:
    -------
    ValueError
        If the `model_source` is not 'GPT' or 'LLAMA'.
    """
    
    if model_source == 'GPT':
        url = "https://api.openai.com/v1/chat/completions"
        model = "gpt-4o-mini"
    elif model_source == 'LLAMA':    
        url = "https://apphubai.wolke.uni-greifswald.de/v1/chat/completions"
        model = 'gpt-3.5-turbo'
    else:
        raise ValueError(f'Variable {model_source=} must either be "GPT" or "LLAMA"')

    # Define the headers
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    

    data = {
        "model": model,
        "messages": [
            {
                "role": "system",
                "content": message_system
            },
            {
                "role": "user",
                "content": message_user
            }
        ]
    }

    # Send the POST request
    response = requests.post(url, headers=headers, json=data)
    llm_response =''
    # Check if the request was successful
    if response.status_code == 200:
        response_data = response.json()
        # Print the response content
        llm_response = response_data['choices'][0]['message']['content']
        print(llm_response)
        print("I delivered!")
    else:
        print(f"Request failed with status code {response.status_code}: {response.text}")

    return llm_response

def main():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description='Query LLM to identify product names that contain plurals or repetitions of words. LLAMA only works with active VPN connection. Note that LLAMA performs very poorly. gpt4o-mini overcalls plurals, but it is feasible to through the output.')

    # Add the input file argument
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Path to the input text file with one product name per line.')

    # Parse the arguments
    args = parser.parse_args()

    # Get the input file path
    input_file = args.input

    # Process the input file (placeholder for the actual functionality)
    with open(input_file, 'r') as file:
        product_names = file.readlines()
        # remove newline characters
        product_names = [name.strip() for name in product_names]
    
    # Particion product names in batches
    nproducts = 5000
    product_names = [', '.join(product_names[i:i+nproducts]) for i in range(0, len(product_names), nproducts)]

    # Placeholder: Code to query LLAMA at the University of Greifswald
    for product in product_names:
        product = product.strip()
        #print("Trying to send query to LLM.")
        response = message_llm(product)
        print(response)
        # ask for pressing enter to continue
        #input("Press Enter to continue...")
        #print("Continuing...")
    print("Finished processing all product names.")


if __name__ == '__main__':
    main()
