
import streamlit as st
from Bio import SeqIO
import io
import tensorflow as tf
import numpy as np

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

st.title('CrisprHAL Online model attempt 2, please work')
uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta","txt"])

if uploaded_file is not None:
    fasta_content = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    targets = find_sgRNA_sequences(fasta_content)
    if targets:
        st.write("Found sgRNA Targets:")
        st.write(targets)
        modelinputs = convert_base(targets)
        #st.write(modelinputs)
        st.write("Running model, please wait...")
        model = load_model()
        predictions = model.predict(modelinputs)
        for i in range(0,len(targets)):
            st.write(str(targets[i]) + "\t" + str(predictions[i]))
    else:
        st.write("No valid sgRNA targets found.")