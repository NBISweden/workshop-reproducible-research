# Tools for Reproducible Research

!!! attention
    If you're viewing this page from the course documentation at Read the Docs
    [document at HackMD](https://hackmd.io/4hSZFZdiROGOq8ZtiAvMUw).
    then you won't be able to add questions/feedback nor see the live updates
    to the document. To do that, go to the

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
- Leif Väremo

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

**Question:**

**Response:**

### Conda

**Question:** I have tried to use conda on Uppmax before, but I haven't been able to get it to work with slurm job files. I always had to activate the environment before submitting the job, rather than doing it within the job. Is there any way around this? This is the error I get when trying to activate conda within a bash script: ![](https://i.imgur.com/zbs9qXU.png)

**Response 1:** Yes, it is possible to activate conda environments within bash scripts. You need to add `conda activate environment` to the beginning of the bash script (basically, replacing the line `module load toolA toolB toolC`, with these tools being in the conda environment). You need to make sure the shebang (first line of the script) such as `#!/bin/bash` does _not_ include `-l`. Alternatively, it might work to add the line `source ~/miniconda3/etc/profile.d/conda.sh` before the line that activates the conda environment - this has solved the same error for me on a different HPC cluster (make sure to update the path to `etc/profile.d/conda.sh` depending on your installation directory).

**Response 2:** While you can activate environments within bash scripts and SLURM jobs, it can sometimes lead to issues (like the one in your screenshot; initialising Conda shells is slightly different from within SLURM jobs and similar login-environments; [YMMV](https://en.wiktionary.org/wiki/your_mileage_may_vary)). Activating an environment before submitting jobs to SLURM means that SLURM will use the environment at the point of entry, which is a feature of SLURM. This means that anything you have access to in terms of software and packages at the time of submission is what's available for the SLURM job.

--

**Question:** I have some questions related to the general usage: Is Conda limited to the terminal, including for the packages that are used, or is it possible to combine with external GUIs like RStudio? If it can, can it talk to for example a Windows installation of RStudio or would it need to be installed within the Linux subsystem when using a Windows-computer? Are there other similar softwares that work on Windows? Are there any GUI implementations of Conda?

**Response 1:** Yes some of the tools we've mentioned have built in support for Conda. PyCharm for instance lets you define and create conda environments for new projects that you start from within PyCharm. Jupyter, which we'll cover later in the course, also has some plugins you can install to manage conda environments (in that case from your browser). 

Regarding Rstudio you can install Rstudio with conda (`conda install -c r rstudio`). We've found that there are sometimes issues with making the conda-installed version of Rstudio aware of other R-related paths on your system though.

**Response 2:** It should be possible to activate a conda environment on the Anaconda prompt on your Windows PC and then start RStudio from there with `call rstudio.exe` (see [here](https://community.rstudio.com/t/opening-r-studio-with-conda-env-version-of-r-fails-unless-launched-via-conda-prompt/72305)). 


### Conda tutorial follow-up questions

_Write your questions about Conda here that come up later during the week._

**Question:** PyCharm looks great but a bit overwhelming, is there any good website for starting with it?

**Response:** There is an introduction at the [PyCharm website](https://www.jetbrains.com/help/pycharm/quick-start-guide.html#ui) that might be suitable. You can also search on YouTube, which has a few videos as well.

## Day 3 (Wed Nov 25)
## Day 4 (Thu Nov 26)
## Day 5 (Thu Nov 27)