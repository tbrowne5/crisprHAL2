
import streamlit as st
from Bio import SeqIO
import io
import tensorflow as tf
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde

def load_model():
    model = tf.keras.models.load_model('Citro_TevSpCas9.h5')
    return model

def convert_base(arrayinput):
    arrayoutput = []
    for sequence in arrayinput:
        onehotencoding = []
        for i in range(len(sequence)):
            if sequence[i].upper() == "A":
                onehotencoding.append([1,0,0,0])
            elif sequence[i].upper() == "C":
                onehotencoding.append([0,1,0,0])
            elif sequence[i].upper() == "G":
                onehotencoding.append([0,0,1,0])
            elif sequence[i].upper() == "T":
                onehotencoding.append([0,0,0,1])
            elif sequence[i].upper() == "N":
                onehotencoding.append([0,0,0,0])
        arrayoutput.append(np.array(onehotencoding))
    return np.array(arrayoutput)

def find_sgRNA_sequences(fasta_content):
    # Assuming the sgRNA targets are to be 20 nucleotides long
    target_length = 28
    targets = []

    for record in SeqIO.parse(fasta_content, "fasta"):
        sequence = str(record.seq)
        # Sliding window to find all possible sgRNA targets of specified length
        for i in range(len(sequence) - target_length + 1):
            sgRNA = sequence[i:i + target_length]
            # Example condition: sgRNA must start with 'G'
            if sgRNA[21:23] == "GG":
                targets.append(sgRNA)
        sequence = str(record.seq.reverse_complement())
        # Sliding window to find all possible sgRNA targets of specified length
        for i in range(len(sequence) - target_length + 1):
            sgRNA = sequence[i:i + target_length]
            # Example condition: sgRNA must start with 'G'
            if sgRNA[21:23] == "GG":
                targets.append(sgRNA)

    return targets

# CURRENTLY TESTING THIS FUNCTION TO DOWNLOAD FILES...
def output_processing(np_train, y_pred):
    np_train = np.array(np_train)
    outputreturn = {}
    for i in range(0,len(np_train)): outputreturn[np_train[i][0:20]] = str(format(float(y_pred[i]), '.8f')).replace(' [','').replace('[', '').replace(']', '')
    return dict(sorted(outputreturn.items(), key=lambda item: float(item[1]), reverse=True))

def output_download(df):
    
    # Sort DataFrame in descending order by the first column as an example
    #sorted_df = df.sort_values(by=0, ascending=False)
    # Convert DataFrame to CSV
    csv = df.to_csv(index=True)
    # Create a download link
    st.download_button(
        label="Download CSV",
        data=csv,
        file_name='sorted_sgRNA_targets.csv',
        mime='text/csv',
    )

def generate_plot(preds):
    data = np.array(preds) #np.random.normal(0, 1, size=100)
    #data = np.sort(data)
    #data = data[::-1]
    st.write(np.max(data))

    # Create a histogram using Matplotlib
    #fig, ax = plt.subplots()
    #ax.hist(data, bins=20, alpha=0.5, color='blue')
    #ax.set_title('Histogram (Matplotlib)')
    #ax.set_xlabel('Values')
    #ax.set_ylabel('Frequency')

    # Display the plot
    #st.pyplot(fig)

    # Create a density plot using Seaborn
    fig2, ax2 = plt.subplots()
    backgroundalpha = 0.3
    ax2.axvspan(-3.5, -2.5, color='#D9DD17', alpha=backgroundalpha)    # Coloring the background from -3 to -2
    ax2.axvspan(-2.5, -1.5, color='#A0D52F', alpha=backgroundalpha) # Coloring the background from -2 to -1
    ax2.axvspan(-1.5, -0.5, color='#70C94B', alpha=backgroundalpha)  # Coloring the background from -1 to 0
    ax2.axvspan(-0.5, 0.5, color='#4ABB5F', alpha=backgroundalpha)    # Coloring the background from 0 to 1
    ax2.axvspan(0.5, 1.5, color='#2DA86F', alpha=backgroundalpha)    # Coloring the background from 0 to 1
    ax2.axvspan(1.5, 2.5, color='#1D9677', alpha=backgroundalpha)     # Coloring the background from 1 to 2
    ax2.axvspan(2.5, 3.5, color='#1B847D', alpha=backgroundalpha)   # Coloring the background from 2 to 3
    sns.histplot(data, kde=True, color='black', ax=ax2)
    labels = ["Rare","Very Low","Low","Moderate","High","Very High","Rare"]
    for i in range(0,7):
        ax2.annotate(labels[i], xy=(((i)/7)+(1/14), 0.915), xycoords='axes fraction', xytext=(0, 10), 
             textcoords='offset points', ha='center', va='bottom',
             fontsize=8, color='black')
    plt.axvline(x=float(np.max(data)), color='red', linestyle='--', linewidth=2)
    plt.annotate('Best', xy=(float(np.max(data))+0.4, 0.9), xycoords='data', xytext=(0, 10), 
             textcoords='offset points', ha='center', va='bottom',
             fontsize=12, color='red')
    ax2.legend([],[], frameon=False)
    ax2.set_xlim(-3.5, 3.5)
    ax2.set_title('How to interpret predicted scores')
    ax2.set_xlabel('Predicted activity scores')
    ax2.set_ylabel('Density')

    # Display the plot
    st.pyplot(fig2)

st.title('crisprHAL: bacterial Cas9/sgRNA activity prediction')
uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta","txt"])

if uploaded_file is not None:
    fasta_content = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    targets = find_sgRNA_sequences(fasta_content)
    if targets:
        #st.write("Found sgRNA Targets:")
        #st.write(targets)
        modelinputs = convert_base(targets)
        #st.write(modelinputs)
        #st.write("Running model, please wait...")
        model = load_model()
        predictions = model.predict(modelinputs)
        #predictionsdf = pd.DataFrame(predictions, index=predictions[:,0])
        output = output_processing(targets,predictions)
        
        st.write("Best guides found:")
        top=1
        for item in output:
            st.write("#" + str(top) + " sgRNA: " + str(item[0:20]) + ", score: " + str(output[item]))
            top+=1
            if top > 5: break
        #st.write("\n\n\n\n\nALL OUTPUTS:")
        #for item in output:
        #    st.write(" sgRNA: " + str(item[0:20]) + ", score: " + str(output[item]))
        outputdf = pd.DataFrame.from_dict(output,orient='index')
        #st.write(outputdf)
        output_download(outputdf)
        generate_plot(predictions)
    else:
        st.write("No valid sgRNA targets found.")
