# Tools for Reproducible Research

!!! attention
    If you're viewing this page from the course documentation at Read the Docs, then you won't be able to add questions/feedback nor see the live updates
    to the document. To edit go to the [document at HackMD](https://hackmd.io/4hSZFZdiROGOq8ZtiAvMUw).

This is a collaborative document where you can give feedback, ask questions and receive help from both TAs and other students. The idea is that it will complement the Zoom breakout rooms and allow you to get more rapid help with common questions and issues.

We've created headers for each day to keep things organised. Simply write your questions and comments under the respective day.

## How to edit

Simply click the `Edit` button in the top menu to get started. You can view the document in split Markdown, Markdown/Preview or Preview modes.

## Important links

- Zoom link: [https://stockholmuniversity.zoom.us/j/62925138459](https://stockholmuniversity.zoom.us/j/62925138459)
- [Course homepage](https://nbis-reproducible-research.readthedocs.io/en/latest/)

## Teachers

- John Sundh
- Verena Kutschera
- Erik Fasterius
- Tomas Larsson
- Leif Väremo Wigge

## Technical issues

**Problem:** Conda crashes when trying to install the snakemake environment from `snakemake/environment.yml`.

**Solution:**  Install snakemake manually by typing `conda create -n snakemake -c bioconda snakemake`. Activate the environment then install `wget` and `graphviz`: `conda install graphviz wget`.

## Day 1 (Mon Nov 23)

### Git

**Question:** During setup, we changed the default name for the main git branch to `main`. In the beginning of the tutorial, two people in our group got the message that they were on branch `master` when typing `git status`, instead of branch `main`. They are currently updating their git installation as both had older versions. Do you have any other idea what might have gone wrong here?

**Response:** This is most likely because of the older Git version: the parameter `init.defaultBranch` that is used was introduced in version 2.28, so any version before that would just have created a `master` branch, even if that parameter was set to `main`.

--

**Syntax for annotating tags**: Since it's not really clear from the tutorial, here's an example of how to annotate a tag upon creation: `git tag -a submission1 -m "Tag for first submission"`. If you don't use the `-m` flag you will open up your default editor and can write the annotation there (exactly in the same way as for commit messages), which is generally recommended.

--

Some parts of the tutorial is mentioning the `git reset` command, while those of you with an updated Git version will likely have `git status` mentioning `git restore` instead. This is also due to version differences, where the old `git reset` did several things that have now been divided up into some new and more streamlined commands, `git restore` being one of them. Here is a link to the official Git documentation on these two commands (and `git revert`): https://git-scm.com/docs/git#_reset_restore_and_revert

## Day 2 (Tue Nov 24)

### Git

**Question:** What is the difference between `git fetch` and `git pull`?

**Response:** Copy-and-pasted from the text in the tutorial: "Another command is `git fetch`, which will download remote changes without merging them. This can be useful when you want to see if there are any remote changes that you may want to merge, without actually doing it, such as in a collaborative setting. In fact, `git pull` in its default mode is just a shorthand for `git fetch` followed by `git merge FETCH_HEAD` (where `FETCH_HEAD` points to the tip of the branch that was just fetched)." You could thus say that `git fetch` means "download changes from the remote" whereas `git pull` is "download *and merge* changes from the remote".

--

**Comment:** As part of Git tutorial, it would be nice to have a working example on how it works in real life. Like what you guys normally do in your regular day. 

**Response:** The tutorial teaches already most of what we do in our daily work. An easy way to get started if you work alone is to track your code with Git using `git add` and `git commit`, and then to push it to your GitHub repository with `git push`. If you collaborate with others, it is highly recommended to start your work on your code with `git pull` (or `git fetch` followed by `git merge`) to update your local copy of the repository with changes your collaborators might have made. 

--

**Question:** Anything to think about when using git from a server *e.g.* Uppmax? What about servers like Bianca that don't have a direct connection to the internet?

**Response:** Git is a very useful tool to keep identical copies of your code on your computer and on Uppmax. You can push your changes made on your local computer to your remote Git repository and pull them to Uppmax (and vice-versa). For Bianca you will have to transfer the repository via the wharf. One way to do it is to start the repository from scratch on Bianca, work with it as you would typically by adding commits etc. Then when a project is finished you transfer the git directory (excluding sensitive data of course) via the wharf to another location, from where you can also push the repository to a remote.

--

**Question:** Are there graphical interfaces for Git?

**Response:** Yes, there is [GitHub Desktop](https://desktop.github.com/), for example. But programs like PyCharm and RStudio usually also have Git so you can add, commit and push your changes from within these programs. If you're one of the people who like Vim there's plugins like `vim-gitgutter` and `vim-fugitive`, which are excellent.

### Git tutorial follow-up questions

_Write your questions about Git here that come up later during the week._

**Question:** Hi, in section the “Remote branches”, I didn’t manage to push changes after merging trimming branch to the main branch. I did:

```bash
git checkout master \
git merge trimming \
git push origin master \
Username for 'https://github.com': <username>
Password for 'https://emiuga@github.com': <password>
To https://github.com/emiuga/git_tutorial.git
 ! [rejected]        master -> master (non-fast-forward)
error: failed to push some refs to 'https://github.com/emiuga/git_tutorial.git'
hint: Updates were rejected because the tip of your current branch is behind
hint: its remote counterpart. Integrate the remote changes (e.g.
hint: 'git pull ...') before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
```

What am I missing? Thanks!

**Response:** The hint is a good place to start here: it tells you that your current branch is behind the remote one, so you have pushed something else without updating the local branch first. Try to do what it tells you, *i.e.* do a `git pull` first (followed by resolving any potential conflicts) and then try pushing again.

--

### Conda

**Question:** I have tried to use conda on Uppmax before, but I haven't been able to get it to work with slurm job files. I always had to activate the environment before submitting the job, rather than doing it within the job. Is there any way around this? This is the error I get when trying to activate conda within a bash script: ![](https://i.imgur.com/zbs9qXU.png)

**Response 1:** Yes, it is possible to activate conda environments within bash scripts. You need to add `conda activate environment` to the beginning of the bash script (basically, replacing the line `module load toolA toolB toolC`, with these tools being in the conda environment). You need to make sure the shebang (first line of the script) such as `#!/bin/bash` does _not_ include `-l`. Alternatively, it might work to add the line `source ~/miniconda3/etc/profile.d/conda.sh` before the line that activates the conda environment - this has solved the same error for me on a different HPC cluster (make sure to update the path to `etc/profile.d/conda.sh` depending on your installation directory).

**Response 2:** While you can activate environments within bash scripts and SLURM jobs, it can sometimes lead to issues (like the one in your screenshot; initialising Conda shells is slightly different from within SLURM jobs and similar login-environments; [YMMV](https://en.wiktionary.org/wiki/your_mileage_may_vary)). Activating an environment before submitting jobs to SLURM means that SLURM will use the environment at the point of entry, which is a feature of SLURM. This means that anything you have access to in terms of software and packages at the time of submission is what's available for the SLURM job.

--

**Question:** I have some questions related to the general usage: Is Conda limited to the terminal, including for the packages that are used, or is it possible to combine with external GUIs like RStudio? If it can, can it talk to for example a Windows installation of RStudio or would it need to be installed within the Linux subsystem when using a Windows-computer? Are there other similar softwares that work on Windows? Are there any GUI implementations of Conda?

**Response 1:** Yes some of the tools we've mentioned have built in support for Conda. PyCharm for instance lets you define and create conda environments for new projects that you start from within PyCharm. Jupyter, which we'll cover later in the course, also has some plugins you can install to manage conda environments (in that case from your browser). 

Regarding Rstudio you can install Rstudio with conda (`conda install -c r rstudio`). We've found that there are sometimes issues with making the conda-installed version of Rstudio aware of other R-related paths on your system though.

**Response 2:** It should be possible to activate a conda environment on the Anaconda prompt on your Windows PC and then start RStudio from there with `call rstudio.exe` (see [here](https://community.rstudio.com/t/opening-r-studio-with-conda-env-version-of-r-fails-unless-launched-via-conda-prompt/72305)).  
_Update on response:_ On the Linux Bash Shell, you can activate a Conda environment (that includes RStudio such as our RMarkdown Conda environment) and start RStudio by simply typing `rstudio &`. This RStudio installation has access to all packages installed in the Conda environment.  
In my case, I had to install and configure a few additional programs to make it work: 1. On the Ubuntu app, a plugin was missing that I found by typing `sudo apt-cache search qt5 | grep 'X11 extras'`, I installed the tool by typing `sudo apt-get install libqt5x11extras5`; 2. I had to install VcXsrc (e.g. see [here, scroll down to "Graphical Applications"](https://seanthegeek.net/234/graphical-linux-applications-bash-ubuntu-windows/) except that I had to add ` export DISPLAY=$(cat /etc/resolv.conf | grep nameserver | awk '{print $2}'):0` to `.bashrc`.

### Conda tutorial follow-up questions

_Write your questions about Conda here that come up later during the week._

**Question:** PyCharm looks great but a bit overwhelming, is there any good website for starting with it?

**Response:** There is an introduction at the [PyCharm website](https://www.jetbrains.com/help/pycharm/quick-start-guide.html#ui) that might be suitable. You can also search on YouTube, which has a few videos as well.

## Day 3 (Wed Nov 25)

### Snakemake

**Question:** I would like to implement these steps for creating an AWS website and uploading images to the website. Any special tips for incorporating AWS into the Snakefile? The PlotCritic directory and scripts all come inside a Singularity container.

```bash
python PlotCritic/setup.py -p NA12878_trio_22 -e ryan@layerlab.org -a YOUR_AWS_ACCESS_KEY_ID -s "YOUR_AWS_SECRET_ACCESS_KEY"

python upload.py -d sv_imgs -c PlotCritic/config.json

python PlotCritic/retrieval.py -c PlotCritic/config.json > retrieved_data.csv

python annotate.py -s retrieved_data.mod.csv -v NA12878.22.vcf -a NA12878.22.score.vcf -o mean -n 1,0,1
```

**Response:** Interesting question. I don't know anything about AWS, but I guess this is all standard except for supplying the AWS access key? I would probably keep the access key on your computer in file with restricted access (only readable/writable by you), then supply the key from the rule parameters, _e.g._: 

```
params:
    AWS_key = config["AWS_key"]
```

Then you add some python code to your Snakefile which reads the key from file and adds it to the config dictionary. But this is just off the top of my head with very little knowledge about AWS :).


--

**Question:** For this rule, could one also change the shell command to the following?

```python
rule concatenate_files:
    input:
        first="{first}.upper.txt",
        second="{second}.upper.txt"
    output:
        "{first}_{second}.txt"
    shell:
        """
        cat {input.first} {input.second} > {output}
        """
```

**Response:** Yes, looking at it it appears to be correct and should work.

--

**Question:** There's a warning from HTSeq about missing index files for the bam files, but the workflow still finishes. Why is this?

**Response:** It appears that HTSeq expects indexed bam files but it's unclear what it does with them (if anything) as at least in this tutorial the output is the same with and without the indices. For future uses of the program it is probably a good idea to supply the index files as well.

--

**Question:**
Just a bit confused about the Logs question, where it says to write the logs to `results/logs/rule_name/`. Should that be a separate `output:` ? I tried throwing it into `2> logs/rule_name/{log}` but that didn't work :P

**Response:** The `log` directive works similar to `input` and `output` in that you specify the path to a file to which the standard error/out (log) should be written to:

```python
rule toolX:
    input: "data/input/{sample}_file"
    output: "results/toolX/{sample}_file_processed.out"
    log: "results/logs/toolX/{sample}.log"
    shell:
        """
        toolX {input} > {output} 2> {log}
        """
```

The `log` directive exists separately from `output` because output files are automatically deleted by snakemake if the rule fails for some reason. The file specified with `log` will be kept even when the rule fails, so that you can find out what the standard error/out of that rule was and can (hopefully) deal with the error.

By specifying the full path to the log file with the `log` directive, snakemake will find the correct log file to keep when the rule fails. In your solution where you specified the path to the log file in the `shell`, snakemake searches for the file you specified with `log` (in my example that would be `"{sample}.log"`) in the main snakemake directory because it would not know the full path.

--

**Question:** Can we have an example that we can go back to reference of how to use snakemake on Uppmax?

**Response:** The relevant parameter you need to use is `--cluster` (read more [here](https://snakemake.readthedocs.io/en/v5.30.1/executing/cli.html?highlight=sbatch#CLUSTER)). You can then use a script to run snakemake in cluster execution mode, which will then send jobs off to Uppmax (SLURM) instead of running them locally. Snakemake will take care of the SLURM queue, monitoring jobs, *etc.*, and there's even a Uppmax Snakemake profile you can use. An example is this:

```bash
snakemake \
    --cluster "sbatch \
                  --account {cluster.account} \
                  --partition {cluster.partition} \
                  --ntasks {cluster.ntasks} \
                  --time {cluster.time}
```

The SLURM variables you normally write in the beginning of your script is thus given at the command line. The above example uses a `cluster.yml` where the values are stored (which works similarly to the `config.yml` file used in the tutorials). Additionally, you also need to use a session manager like [tmux](https://github.com/tmux/tmux/wiki), so that you don't get logged out of Uppmax while Snakemake is running.

**SLURM Profile**: Another option is the [SLURM profile](https://github.com/Snakemake-Profiles/slurm) for Snakemake. The way you use it is by 'deploying' it using the software `cookiecutter` (install through conda), at setup it asks you for the project account id to use for submitting jobs. See the [README](https://github.com/Snakemake-Profiles/slurm#example-1-project-setup-to-use-specific-slurm-account) for how to set up with a specific slurm account.

Then you supply information on how to submit rules via the `resources:` directive, _e.g._:

```python
rule align_to_genome:
    input:
        "{genome_id}.bt2",
        "{sample}.fastq.gz"
    output:
        "{sample}.bam"
    resources:
        runtime = 360
    threads: 10
    shell:
        """
        aligner -t {threads} -i {input[1]} -x {input[0]} > {output}
        """
```

Then run the workflow as:
```bash
snakemake -j 10 --profile slurm
```

Any rule with `runtime` in the `resources` directive will be submitted to the queue with `runtime` as the allocated time, and `{threads}` as the allocated cores.

You can even specify if rules should be resubmitted to the queue asking for more resources on subsequent attempts. Do this by modifying the `resources` directive with _e.g._:

```python
    resources:
        runtime = lambda wildcards, attempt: attempt*360
```

This way the rule will ask for 360 minutes the first attempt, 2*360 minutes on the second attempt etc.

Also see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#standard-resources) about specifying resources.

--

**Question:** How do you run Snakemake with singularity? Follow-up: is it possible to specify the path to a `*.sif` file with the `singularity` directive?

**Response:** You can use the `singularity` directive equivalently to the `conda` directive in your snakemake rule to specify a Singularity or Docker image:

```python
rule sort_bam:
    """
    Sort a bam file.
    """
    input:
        "intermediate/{sra_id}.bam"
    output:
        "intermediate/{sra_id}.sorted.bam"
    singularity:  "docker://biocontainers/samtools:v1.9-4-deb_cv1"     
    shell:
        """
        samtools sort {input} > {output}
        """
```

When you execute the workflow, you will do that by providing a flag that tells snakemake to look for the `singularity` directive and to pull the image specified with the directive from e.g. DockerHub:

```bash
snakemake --use-singularity
```

**Response to follow-up**: Yes! It is also possible to use a `*.sif` file by providing the path directly in the `singularity: ` directive, _e.g._:

```python
    singularity: "path/to/singularityfile.sif"
```


--

### Snakemake tutorial follow-up questions

_Write your questions about Snakemake here that come up later during the week._

**Question:**

**Response:**

--

## Day 4 (Thu Nov 26)

### Jupyter

**Question:** How do you launch a Jupyter notebook from UPPMAX? Is there a module you can load?

**Response 1:** You can try the [ThinLinc](https://www.uppmax.uu.se/support/user-guides/thinlinc-graphical-connection-guide/) software that's installed on Uppmax, which can give you a GUI desktop on a login-node on both Rackham and Bianca. None of us teachers have used it, so we can't say any details on it, but we think it might be able to run Jupyter on Uppmax. You can always ask the Uppmax Support for help regarding this, if you're interested in getting it to work.

**Response 2:** On Rackham (but probably not Bianca) you can follow the instructions from the tutorial:

```bash
ssh me@rackham.uppmax.uu.se -L8888:localhost:8888
jupyter notebook --ip 0.0.0.0 --no-browser
```

--

**Question:** What does the icon in the upper right corner mean that says "Not trusted" or "Trusted"?

**Response:** The notebook gets "Trusted" when the user has executed all code, and whenever you load it again, HTML and Javascript output will be displayed. For "Not trusted" notebooks, these are disabled until you run the cells that require either of them. You can read more about it here: https://jupyter-notebook.readthedocs.io/en/latest/notebook.html#trusting-notebooks

--

**Question:** I tried the widget "ColorPicker", but the color won't change from "red" to "blue". Why could that be? 

```python
colorpicker = widgets.ColorPicker(
    description='Color',
    value='blue')
```

**Response:** You also need to edit the following line of code from:

```python
interactive_plot = interactive(sine_curve, A=A, f=5, p=5)
```

to

```python
interactive_plot = interactive(sine_curve, 
            A=(1, 5, 1), 
            f=(0, 5, 1), 
            p=(1, 5, 0.5),
            color=colorpicker) ## <- Supply the colorpicker to the function
```

### R Markdown

**Question:** Would it be possible to give an example on how you start up a new R Markdown document on uppmax, using your preferred editor (whatever that might be)?

**Response:** As R Markdown is simply a format, you can just create a new file ending in `.Rmd` and start typing away. Depending on your editor this can require some sort of port/X-forwarding or GUI interface with [ThinLinc](https://www.uppmax.uu.se/support/user-guides/thinlinc-graphical-connection-guide/) or similar, but it can also not matter if you're using a text editor like Vim or Emacs.

With Vim, for example, I (Erik) has so-called "snippets" that allow me to easily start from a fresh R Markdown file and quickly add a template YAML header (which I edit) and then start typing away at markdown and chunks as needed. I render with the same command line-based code that's shown in the tutorial (through Vim shortcuts). If you're interested I could show in more detail, but it varies very much from setup to setup. For RStudio it's the same as editing any other R Markdown document, you just have to connect it to Rackham through something. This can be done if you have X-forwarding enabled when you SSH in (read more at the [Uppmax website](https://www.uppmax.uu.se/support/user-guides/r-user-guide/)).

A thing to keep in mind with ANY notebook that you run on Uppmax (so either Jupyter or R Markdown) is that the login nodes are not capable of very heavy computations (probably less than your own laptop, actually), so if you're doing heavy computations you'll want to submit it to the SLURM system (either through scripts or, more reproducibly, through a workflow manager like Snakemake).

### Docker

**Question:** Do images you build become independent from the base image layer you use in the `FROM` directive?

**Response:** Yes, once an image is built it is a separate stand-alone unit of software. At build-time however the base image you specify in `FROM` has to be available to Docker somehow, either via a remote registry such as DockerHub or on your local computer.

--

**Question:** Can you push images that you haven't built with a tag (`-t`) to a remote repository?

**Response:** Yes but you have to tag the image first with `docker tag`.

### Singularity

**Comments for Windows users:** 
- The website has been updated with instructions how to open the Vagrant VirtualBox (https://nbis-reproducible-research.readthedocs.io/en/latest/singularity/#setup).
- Files in the Vagrant VirtualBox and your local computer are shared via the folder `vagrant/` (in the Vagrant VirtualBox). Per default, you are standing in `/home/vagrant`.
 
**Comment for Linux users** (Ubuntu 16.04)
I run into trouble with version 2.6.0 as described in the setup. I had to install the newest by following [this link](https://sylabs.io/guides/3.4/user-guide/installation.html#distribution-packages-of-singularity) (but actually installing 3.7.0 as it is the newest).  Make sure to install Go as described before singularity.

### Jupyter tutorial follow-up questions

_Write your questions about Jupyter here that come up later during the week._

**Question:**
I had to create a separate conda environment for jupyter notebook on Uppmax:

```yaml
channels:
- conda-forge

dependencies:
- jupyter
- nb_conda
- matplotlib
- ipywidgets
- seaborn
- pandas
```

Then I activate this: `conda activate jupyter`

But when I run `jupyter notebook --ip 0.0.0.0 --no-browser`, the link provided doesn't seem to work :)

   To access the notebook, open this file in a browser:
        file:///domus/h1/.local/share/jupyter/runtime/nbserver-36818-open.html

**Response:** We'll look into that :-) 

--

### R Markdown tutorial follow-up questions

_Write your questions about R Markdown here that come up later during the week._

**Comment:** A comment on writing simple reports: If you write scripts in R and use RStudio, you can also simply knit R script without using R markdown language, as RStudio converts a script for you.

--

### Docker tutorial follow-up questions

_Write your questions about Docker here that come up later during the week._

**Question:**

**Response:**

--

### Singularity tutorial follow-up questions

_Write your questions about Singularity here that come up later during the week._

**Question:**

**Response:**

--

## Day 5 (Thu Nov 27)

### Docker

**Comment:** Some Mac users might run into an error message when building the Docker image for the MRSA case study from `Dockerfile`, during the step when conda is installed. This can be solved by giving the Docker Desktop App more memory, as it only has access to 2 GB RAM per default. Go to the Docker Desktop App, click on the little wheel in the upper right corner, and go to "Resources". Under "Advanced", increase the memory to *e.g.* 4 GB and click "Apply & Restart".

### Singularity

--

### Putting the pieces together

_Write your questions about the lecture here_

### Q&A: How to implement these procedures on a day-to-day basis

**Question:** How do you deal with files not neccesarily following the same naming convention? For example, I work with genome assemblies from different people and different tools, and the headers are often very different in regards to which chromosome a scaffold represents, so I need different code for each assembly to change the header to useful information for my finaly results.

**Response:** You are sadly correct that if you have input data that is used for the same thing (*e.g.* alignment) but are formatted differently, you need to have some kind of variation on those files are handled. You could use separate Snakemake rules that take the differently formatted+named files, for example, or you could have a single rule with some script with if-clauses.

--

**Question:** Will we get a certificate for this course?

**Response:** Yes!

--

**Question:** Mamba seems good too? Any thoughts?

**Response:** It's great, use it!

--

**Question:** Any suggestion about a user-friendly text editor in Ubuntu? e.g. that would make easier to work more interactively and documenting at the same time, similar to Rstudio.

**Response:** How about Rstudio? There's an Ubuntu version of that. Sublime, Atom, VSCode, PyCharm, Visual Studio Code and Vim are others.

--

**Question:** Technical question: do you need to supply a .sif file in this command?

```bash
snakemake -s snakefile --use-singularity *.sif -j 1
```

or just

```bash
snakemake -s snakefile --use-singularity -j 1
```

**Response:** The `--use-singularity` flag only tells snakemake to use singularity for the rules that have the `singularity` (more recently `container`) directive specified. Or for workflows that have a global `singularity` directive. So the latter command is correct. A rule using singularity may look like this:

**Note:** Recently snakemake has made the `singularity` directive deprecated. Instead you should use `container`, but this may depend on the snakemake version you're using.

```python
rule trimrule:
    input: "file1.fastq.gz"
    output: "file1.trimmed.fastq.gz"
    singularity: "docker://biocontainers/seqtk"
    shell:
        """
        seqtk trimfq {input} > {output}
        """
```

--

**Question**: If you want to work with bioinformatics at e.g. NBIS, how do you go about to acquire the skills needed (except for practicing a lot)?

**Reponse:** Practice a lot! Seriously. More useful ideas are: take courses (like you're doing now), google a lot, find best practices, ask colleagues, come to the NBIS drop-in, try to add new things to learn for each new project.

--

**Question**: Any suggestions for learning Python? There's only one course per year at NBIS.

**Response 1:** Have a look at the tutorials at [RealPython](https://realpython.com/).

**Response 2:** I heard positive things about the Python course at [kaggle](https://www.kaggle.com/learn/overview). If you are looking for some examples to apply Python too, check the [Rosalind problems](http://rosalind.info/problems/list-view/).

**Response 3:** Some useful things here: (https://www.codecademy.com)


--

**Question**: Do some of you do even stats and plotting in python?

**Response 1:** Yes! I (John) typically do all my plotting and most statistics through python (via Jupyter notebooks). For plotting I really recommend the `seaborn` package. For stats you have `scikit-learn`, `scipy` etc.

**Response 2:** I'm also doing plotting in Python (Verena) and I'm planning to move statistical analysis to Python, too. I'm currently working through this book: https://github.com/wesm/pydata-book

**Response 3:** (from the zoom chat :-) ) While starting my bioinformatic ambition I tried the capstone project with John Hopkins course on coursera, and this course helped me a lot. So if one has some idea about python and start with something as a reference point, one can start with this. https://www.coursera.org/specializations/genomic-data-science

--

**Question**: This has been a steep one for me, because there is ggplot analogues in Python but I already have easy ggplot in R :P will need to spend lots of time I guess! Also seems like more examples for ggplot in R versus nice plotting in Python? Porbably I need to google more :)

**Response:** There is actually a ggplot implementation in Python, but it sadly is not as feature-rich as the one in R and it can be a little buggy. Python plotting libraries are, for example, `bokeh` and `seaborn`.

**Question**: Matplotlib is perhaps less desirable? Thanks again for the great course!

**Response:** Matplotlib is a very powerful but also a very dense package with lots of features and functionality. This can make it somewhat difficult to start using for your daily plotting. Seaborn, for instance, is built on top of Matplotlib but offers a simpler interface to start creating high-quality plots.

--

**Question:** There seems to be an issue with pulling a singularity container through docker, where one needs to specific the version?

For example:

```bash
singularity pull docker://quay.io/biocontainers/*:latest *_latest.sif
```

Or how does one need to specify a specific version?

For example, I have:
docker pull quay.io/biocontainers/smoove:<tag> for tag, can I type: docker pull quay.io/biocontainers/smoove:0.2.6-0 

**Response 1:** I think this might be because a specific image does not have the `latest` tag? Could you specify what image you tried this with? Also, for this example the command should be:

```bash
# For docker
docker pull quay.io/biocontainers/smoove:0.2.6--0
# For singularity
singularity pull docker://quay.io/biocontainers/smoove:0.2.6--0
```

**Response 2:** If you are looking for containers related to bioinformatics, you can also search here: https://biocontainers.pro/#/registry Here, all conda packages are listed in their conda version and their equivalent container version.

--

**Question:** A practical question regarding implimenting these things in our projects. If I want to start implimenting git with a sub-project of what will become a paper (i.e. there are other directories outside this one, but part of the same project), should I make a git for the bigger project now, or can I add this sub-project git to the main project git if I make that at a later date?

**Response:** You can break out part of the project and create a git repository for that (with a paper directory inside it), but try to organise the whole project at one point.