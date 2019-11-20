# Using the tools to put everything together

It is time to try to setup a project from scratch and use the different 
tools that we have covered during the course together! 

This exercise if very open-ended and you have free hands to try out a 
bit what you want. But you should aim to use what you've learned to 
do the following:

1. Create a new git repository for the project (either on BitBucket or GitHub)
2. Add a README file which should contain the required information on how to run the project
3. Create a conda `environment.yml` file with the required dependencies
4. Create a snakemake file to run your workflow
5. Use a `config.yml` file to add settings to the workflow
6. Use git to commit changes for the repository

## Option 1
One option is to try to implement these methods on 
one of your current projects. It is up to you what tools to include in 
making your project reproducible, but aim for at least including git and conda.

!!! tip
    If your analysis project contains 
    computationally intense steps it may be good to scale them down for 
    the sake of the exercise. 

## Option 2
If you don't want to use a project you're currently working on we have 
a suggestion for a small-scale project for you.

The idea is to analyze student experience for this Reproducible Research
course. For this you will use responses from students to the registration 
form for the course. Below you'll find links to **csv** format files
with answers from 3 course instances:

* 2018-11:https://docs.google.com/spreadsheets/d/1yLcJL-rIAO51wWCPrAdSqZvCJswTqTSt4cFFe_eTjlQ/export?format=csv
* 2019-05:https://docs.google.com/spreadsheets/d/1mBp857raqQk32xGnQHd6Ys8oZALgf6KaFehfdwqM53s/export?format=csv
* 2019-11:https://docs.google.com/spreadsheets/d/1aLGpS9WKvmYRnsdmvvgX_4j9hyjzJdJCkkQdqWq-uvw/export?format=csv

The goal here is to create a snakemake workflow which:

1. downloads the csv files (making use of a `config.yml` file to pass the urls)
2. cleans the files (using `wildcards`)

The final step is to plot the student experience in some way.

The first two steps should be part of the workflow. If you need some help
with the cleaning step, see below for a script that you can save to a file
and run on your computer.

??? note "Click to show a script for cleaning column names"
    ```python
    #!/usr/bin/env python
    import pandas as pd
    from argparse import ArgumentParser
    
    def main(args):
        df = pd.read_csv(args.input, header=0)
        df.rename(columns=lambda x: x.split("[")[-1].rstrip("]"), inplace=True)
        df.rename(columns={'R Markdown': 'RMarkdown'}, inplace=True)
        df.to_csv(args.output, index=False)
        
    if __name__ == '__main__':
        parser = ArgumentParser()
        parser.add_argument("input", type=str,
                            help="Input csv file")
        parser.add_argument("output", type=str,
                            help="Output csv file cleaned")
        args = parser.parse_args()
        main(args)
    ```
    
The last step is really up to you how to implement. You could:

* include the plotting in the workflow using an RMarkdown document that 
gets rendered into a report
* have a script that produces separate figures (e.g. `png` files)
* create a jupyter notebook that reads the cleaned output from the workflow
and generates some plot or does other additional analyses

If you need some help/inspiration with plotting the results, click below
to see an example python script that you can save to file and run with
the cleaned files as input.

??? note "Click to show a script for plotting the student experience"
    ```python
    #!/usr/bin/env python

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    mpl.use('agg')
    import pandas as pd
    import seaborn as sns
    import numpy as np
    from argparse import ArgumentParser

    def read_files(files):
        """Reads experience counts and concatenates into one dataframe"""
        df = pd.DataFrame()
        for i, f in enumerate(files):
            # Extract date
            d = f.split(".")[0]
            _df = pd.read_csv(f, sep=",", header=0)
            # Assign date
            _df = _df.assign(Date=pd.Series([d]*len(_df), index=_df.index))
            if i==0:
                df = _df.copy()
            else:
                df = pd.concat([df,_df], sort=True)
        return df.reset_index().drop("index",axis=1).fillna(0)

    def count_experience(df, normalize=False):
        """Generates long format dataframe of counts"""
        df_l = pd.DataFrame()
        for software in df.columns:
            if software=="Date":
                continue
            # Groupby software and count
            _df = df.groupby(["Date",software]).count().iloc[:,0].reset_index()
            _df.columns = ["Date","Experience","Count"]
            _df = _df.assign(Software=pd.Series([software]*len(_df), index=_df.index))
            if normalize:
                _df = pd.merge(_df.groupby("Date").sum().rename(columns={'Count':'Tot'}),_df, left_index=True, right_on="Date")
                _df.Count = _df.Count.div(_df.Tot)*100
                _df.rename(columns={'Count': '%'}, inplace=True)
            df_l = pd.concat([df_l, _df], sort=True)
        df_l.loc[df_l.Experience==0,"Experience"] = np.nan
        return df_l


    def plot_catplot(df, outdir, figname, y, palette="Blues"):
        """Plot barplots of user experience per software"""
        ax = sns.catplot(data=df, x="Date", col="Software", col_wrap=3, y=y, hue="Experience", height=2.8,
                         kind="bar",
                         hue_order=["Never heard of it", "Heard of it but haven't used it", "Tried it once or twice",
                                    "Use it"],
                         col_order=["Conda", "Git", "Snakemake", "Jupyter", "RMarkdown", "Docker", "Singularity"],
                         palette=palette)
        ax.set_titles("{col_name}")
        plt.savefig("{}/{}".format(outdir, figname), bbox_to_inches="tight", dpi=300)
        plt.close()

    def plot_barplot(df, outdir, figname, x):
        """Plot a barplot summarizing user experience over all software"""
        ax = sns.barplot(data=df, hue="Date", y="Experience", x=x, errwidth=.5,
                    order=["Never heard of it", "Heard of it but haven't used it", "Tried it once or twice", "Use it"])
        plt.savefig("{}/{}".format(outdir, figname), bbox_inches="tight", dpi=300)
        plt.close()

    def main(args):
        # Read all csv files
        df = read_files(args.files)
        # Count experience
        df_l = count_experience(df)
        # Count and normalize experience
        df_lp = count_experience(df, normalize=True)
        # Plot catplot of student experience
        plot_catplot(df_l, args.outdir, "exp_counts.png", y="Count")
        # Plot catplot of student experience in %
        plot_catplot(df_lp, args.outdir, "exp_percent.png", y="%", palette="Reds")
        # Plot barplot of experience
        plot_barplot(df_lp, args.outdir, "exp_barplot.png", x="%")

    if __name__ == '__main__':
        parser = ArgumentParser()
        parser.add_argument("files", nargs="+",
                            help="CSV files with student experience to produce plots for")
        parser.add_argument("--outdir", type=str, default=".",
                            help="Output directory for plots (defaults to current directory)")
        args = parser.parse_args()
        main(args)
    ```
    
!!! attention
    Remember to:
    
    1. keep everything versioned controlled with `git` 
    2. add information to the `README` file so others know how to run the project
    3. add required software to the conda `environment.yml` file 