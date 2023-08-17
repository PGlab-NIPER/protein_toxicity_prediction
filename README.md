# protein_toxicity_prediction

Welcome to the Protein Toxicity Prediction Repository.
This prediction script empowers users to assess the potential toxicity of proteins in terms of Cardiotoxicity, Neurotoxicity, and Enterotoxicty. By providing the primary protein sequence of your query protein, our well-trained models will swiftly analyze and predict whether the entered sequence exhibits toxicity, and if so, which specific toxic category it falls under.

# Features of this ML classification model
- Predicts protein toxicity based on Cardiotoxicity, Neurotoxicity, and Enterotoxicty.
- Utilizes advanced machine learning models for accurate predictions.
- Provides results in a convenient .csv format.

# Contents
- Prediction model script written in the R programming language.
- To save you time and resources, we have included the best-performing pre-trained models. These models are saved as pickled files, ensuring you don't need to perform recurrent model training.
- Datasets that contains best features which has highest highimportance in building all three models.
- Datasets that contain the reduced list of features used to build the prediction models after performing feature selection.
  These features have demonstrated the highest importance in creating accurate predictions for all three toxicity categories.
 
# How to Use
1. Make sure you have **R 4.3.0** and **R Studio 2023.06.0** installed on your system.
2. Clone this repository and navigate to the project directory.
3. Open the prediction script in R Studio.
4. Input your protein's primary sequence in the designated area.
5. Run the script to generate predictions.
6. The results will be saved as a .csv file for your reference.
