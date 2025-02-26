# A Bayesian Model of Galactic HII Region Physio-kinematic Distances

## Project Objectives

### 0. Logging research progress

   **Status: Ongoing**
   
   *Learning objective*: Recognize the importance of tracking research objectives, progress, and questions.

   *Criteria for success*: Keep organized, meticulous notes of your research.

   <details>
   The most important part of the research process is probably being able to effectively communicate about the project. This means being able to explain to a random stranger on the street what you're doing, why it's important, and what it means. This is only possible if YOU know what you're doing. To this end, I ask that you keep diligent notes about everything you do related to this project. These notes don't have to be in any specific format, although it would be useful if they were saved in some what that I could also access them (like a google doc). Keep a record of what you do (e.g., I read this paper, I wrote a program that does this, I got confused about this topic, etc.), keep a record of what you want to do next (e.g., I need to write a program that does this other thing, I need to read about this topic, etc.), and, most importantly, keep track of all of the questions that come up (what does this acronym mean, how does this physical thing relate to this other physical thing, etc.). These notes will be invaluable to you as you work on the project. I often get distracted by other tasks and come back to a project after a few days or weeks only to have forgotten what exactly I was doing and what I needed to do next. Without these notes, I would have been lost!
   </details>

### 1. Start writing the "paper"

   **Status: Ongoing**

   *Learning objective*: Develop effective written communication skills

   *Criteria for success*: Start a draft of the introduction/background section of a "paper"

   <details>
   I hope that this project will ultimately result in a publication, but no matter what it will benefit YOU to start writing a "paper" or "final report" for the project right now, before you do anything else. In particular, I want you to focus on the "introduction" section of a paper, where you outline the major research questions and goals of the project. This will immensely benefit you because it will be something that you can look back on when you're knee-deep in data analysis and programming and you've forgotten what the "big picture" of the research project is. Don't worry about the formatting, the specific content, or anything like that now. Just write a paragraph or two about the project, and go back and read/edit it once in a while as you develop a stronger grasp on our research objectives. And it's OK if you don't know what the research questions/goals are yet - that's something we can talk about, which will guide your writing!
   </details>

### 2. Background research

   **Status: Ongoing**

   *Learning objective*: Develop a basic physical understanding of HII regions and the data sets that we'll be working with.

   *Criteria for success*: Be able to draw a cartoon picture of an HII region, relate the physical characteristics of the nebula to its observable properties, and understand how we identify and characterize them.

   <details>
   The first step for any project is to understand what's been done before. In this case, physicists decades ago figured out the basic physics of HII regions and derived all of the equations and relationships that we'll need for this project. Your first challenge is to develop a basic understanding of this work. Here are some resources to get you started, although I hope you will do your own internet-searches to fill in the gaps and answer some questions. Take note of any questions or confusing topics that you come across along the way, and we can talk about them together.

   * Wikipedia: https://en.wikipedia.org/wiki/H_II_region
   * [Essential Radio Astronomy (ERA) Textbook](https://www.cv.nrao.edu/~sransom/web/xxx.html): In particular, chapters 4.2, 4.3, and 7.2 will be useful!
   * [The WISE Catalog of Galactic HII Regions](https://ui.adsabs.harvard.edu/abs/2014ApJS..212....1A/abstract): This paper provides an overview of the latest catalog of Milky Way HII regions. 
   * [The Southern HII Region Discovery Survey](https://ui.adsabs.harvard.edu/abs/2019ApJS..240...24W/abstract): This paper discusses one of many radio recombination line surveys of HII regions (led by yours truly!) 
   </details>

### 3. Preparing your research environment.

   **Status: COMPLETE**

   *Learning objective*: Prepare software environment

   *Criteria for success*: Write a "hello world" program in python

   <details>
   We're going to have to write some computer programs. The first step in this journey will be installing the necessary software on your computer and writing your first python program. There are many ways to set up a python environment, the specifics of which depend on what kind of computer you have, what operating system you use, etc. In general, Google/ChatGPT will probably be more helpful than I. Look up some tutorials, watch some youtube videos, and try to write a "hello world" program in python.
   </details>

### 4. Forward \& Inverse Modeling

   **Status: COMPLETE**

   *Learning objective*: Explore the differences between foward modeling and inverse modeling

   *Criteria for success*: Write a program that can determine the physical conditions of an HII region from observed quantities

   <details>
   Now that you have a "forward model" (a model that can take physical quantities and predict observations), it's time to write an "inverse model" (a model that takes observations and predicts physical quantities). Your objective is to write a program that takes as input the observable properties of an HII region, namely the continuum brightness temperature, the peak radio recombination line brightness temperature, and the FWHM line width, and returns the physical conditions (temperature, emission measure, etc.). Think about how you can adapt the program you've already written in order to achieve this goal! Furthermore, use your existing program to *test* your new program: generate a "synthetic" observation using your existing program and some chosen physical parameters, and then see if your new program is able to predict those physical parameters.
   </details>

### 5. Model Fitting

   **Status: COMPLETE**

   *Learning objective*: Develop an understanding of how models are fit to data

   *Criteria for success*: Implement a least squares algorithm to fit a physical model to some synthetic data

   <details>
   There are a few steps here:
   1. First, you'll need to modify your existing forward model to include *noise* in your simulated spectra. Check out this function for generating "gaussian" noise:
   https://numpy.org/doc/stable/reference/random/generated/numpy.random.normal.html
   2. Next, you'll need to learn a little about the least squares algorithm. Wikipedia is a good start:
   https://en.wikipedia.org/wiki/Least_squares
   3. The most important part of least squares is quantifying the "loss" function, or the "residuals" function. You'll need to write a function that takes as input the model parameters (emission measure, temperature, etc.) as well as the data and then returns the "loss" or the sum of the squared residuals (data - model)^2.
   4. Use your forward model program to generate a synthetic, noisy spectrum.
   5. Finally, use scipy's least_squares implementation to fit your model to the data.
   https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
   </details>

### 6. Monte Carlo Uncertainty Analysis

   **Status: COMPLETE**

   *Learning objective*: Introduce yourself to Monte Carlo methods

   *Criteria for success*: Determine model parameter uncertainties, plot parameter posterior distributions

   <details>
   Monte Carlo analyses are powerful tools for characterizing the uncertainties in model parameters. Use Monte Carlo methods to determine the parameter uncertainties for your forward model. The steps include:
   1. Resampling the data: at each Monte Carlo iteration, resample the data within the assumed uncertainties. This is as simple as adding different random "noise" at each iteration.
   2. Fitting the model: use your least squares implementation to fit the model to the resample data
   3. Parameter posterior: save the best fit parameters at each iteration, then construct histograms of those parameters (the posterior distributions) to characterize the uncertainties
   4. Visualization: plot the posterior distributions, and plot the range of allowed models on the data
   </details>

### 7. Monte Carlo Markov Chain Analysis

   **Status: COMPLETE**

   *Learning objective*: Introduce yourself to Monte Carlo Markov Chain analyses

   *Criteria for success*: Understand Trey's MCMC implementation

   <details>
   Monte Carlo Markov Chain (MCMC) techniques build upon all of the techniques that you've learned so far. Rather than Monte Carlo resampling the data, we now Monte Carlo resample the model parameters from some assumed "prior" distributions, which characterize our prior knowledge of the parameters. Using Markov Chains, we explore the model parameter space in order to construct the "posterior" distributions, the knowledge of the parameters supported by the data. Here are some tasks:
   1. Read about Bayesian statistics: https://en.wikipedia.org/wiki/Bayesian_statistics
   2. Read about MCMC: https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo
   3. Read about `pymc`: https://www.pymc.io/welcome.html
   4. Study the new program: `bayesian_model.py`
   </details>

### 8. Working with Real Data

   **Status: CANCEL**

   *Learning objective*: Understand the data we'll be working with

   *Criteria for success*: Characterize the data, generate summary plots

   <details>
   Finally, we must start exploring the data. You can access the database here: https://doi.org/10.7910/DVN/NQVFLE
   1. Learn about SQL: https://gist.github.com/lionelbarrow/8177236
   2. Learn how to interact with an SQL database in python: https://docs.python.org/3/library/sqlite3.html
   3. Read in the database, and generate plots of HII region properties (positions, LSR velocity, longitude-velocity diagram, etc.)
   </details>

### 9. Testing `physiokinematic`

    **Status: COMPLETE**

    *Learning objective*: Understand the `physiokinematic` forward model

    *Criteria for success*: Run the "test.ipynb" example notebook

    <details>
    I have updated the `physiokineamtic` package and put it on Github: https://github.com/tvwenger/physiokinematic

    Please fork the repository and start investigating the code! Some things to consider:
    1. Try to understand what the programs do
    2. Try to install and run the programs
    3. Specifically, try to run the "test.ipynb" notebook and understand what it is doing
    
    As you work through these steps, note any questions or possible issues that you notice or encounter.
    </details>

### 10. Running `physiokinematic on real data

    **Status: New**

    *Learning objective*: Learn how to process big data sets

    *Criteria for success*: Run the model on all of the data!

    <details>
    Now that we have a working model, it's time to run it on some data!
    Following the process in the "real data" notebook, try fitting the
    model to all of the real HII region data on scylla!

    Take it one step at a time:
    1. Run the model on a single HII region
    2. Save the output (i.e. posterior samples, the "trace") from the model
    3. Create a python script that does (1) and (2)
    4. Update the python script to loop over all HII regions
    5. Compile the data and present the results
    </details>

### 11. Outline Paper

    **Status: New**

    *Learning objective*: Practice writing a scientific manuscript

    *Criteria for success*: Generate a draft outline for the paper

    <details>
    I have shared with you the paper on Overleaf! Currently it
    is filled with a template and demonstration of my writing
    style. Edit this manuscript and create a draft outline for
    our paper! Feel free to follow the template as closely or
    loosely as you like. Also feel free to suggest some catchy
    titles!
    </details>
