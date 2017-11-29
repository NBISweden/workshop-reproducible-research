# git tutorial

## Introduction

Git is a widely used system (both in academia and industry) for version controlling files and
collaborating on code. It is used to track changes in (text) files, thereby establishing a history of all edits made to each file, together with short messages about each change and information about who made it. Git is mainly run from the command line, but there exists several software that have implemented a graphical user interface to run git commands.

Using git to version control and track your files, and edits to those, is an essential step of making your computational research reproducible. A typical git workflow consists of:
* making distinct and related edits to one or several files
* committing those changes (i.e. telling git to add those edits to the history, together with a message about what those changes involve)
* and pushing the commit to a remote repository (i.e. syncing your local project directory with one in the cloud)

There are several benefits with using git in your research project:
* You are automatically forced into a more organized way of working, which is usually a first step towards reproducibility.
* If you have made some changes to e.g. a script and suddenly realize that those were probably not a good idea after all, it is simple to view exactly what changes you have made and revert them.
* If there is more than one person involved in the project, git makes it easy to collaborate on the coding, tracking all edits made by each person, and handle potential conflicting edits.
* Using a cloud based repository hosting service (the one you push your commits to), like e.g. [Github](https://github.com/) or [Bitbucket](https://bitbucket.org/), adds additional features, like being able to comment and discuss the project, code, and edits.
* At some point your project will be published. Github or Bitbucket (or similar) are excellent places to publicly distribute your code. Other researchers can then use git to access your code needed for reproducing your results, in the state it was when used for the publication.
* If needed, you can host private repositories on GitHub and BitBucket as well, which may be convenient during an ongoing research project, before it is publicly published.

The best way to get an idea about what git is simply to start using it. The tutorial below will guide you through the essential steps, with a focus on what is needed for making a project reproducible. There are additional features of both git and the web based repository hosting services (like Github and Bitbucket) that are intentionally not included here. If you are interested of learning more, the web is filled with information (see some examples below)!


#### Tell me more

* For a more complete introduction to git, check e.g.:
https://en.wikipedia.org/wiki/Git
* [A simple git guide](http://rogerdudler.github.io/git-guide/)
* ["Got 15 minutes and want to learn Git?"](https://try.github.io/levels/1/challenges/1)
* [Git reference manual](https://book.git-scm.com/docs)

## Practical exercise

### Setup a new git repository

1. If you have not done so already, go to [bitbucket.org](https://bitbucket.org/account/signup/) and create an account.
* Login and press the plus button to the left and select *Create a new repository*:
  * Make sure you are listed as the owner
  * Add a repository name, e.g. *git_tutorial*
  * You can keep the repo private or make it public, as you wish
  * Skip including a README
  * Make sure Git is selected for version control  
You will now be redirected to the repository page. It is an empty repository, so there is not much to see yet.
* We want to add some content (files) to the repository. To do that we will first *clone* the repository locally:
  * Open a terminal and `cd` to a directory where you want to clone your newly created git repository (perhaps make a new directory for this course). **Important: it should *not* be within the `reproducible_research_course` directory.**
  * Once you are in your directory of choice, run the following command (just make sure to change `user` to your bitbucket username and `git_tutorial` to your repository name, in case you chose something different):
  ```
  git clone git@bitbucket.org:user/git_tutorial.git
  ```
  What will happen now is that the git repository will be cloned (i.e. downloaded) to your computer. You might get a warning that the repository is empty (which in fact is the case).
* A new directory, `git_tutorial` (or a different name if you chose so), has now been created, `cd` into that.  
This is a git version tracked directory. How can you know? Run `git status`! It will probably return something like *On branch master. Initial commit. nothing to commit (create/copy files and use "git add" to track)*. (If you try to run `git status` in a non-git directory, it will say that it is *not a git repository*.) The way this works is that git adds a hidden directory `.git/` in the root of a git tracked directory (run `ls -a` to see it). This hidden directory contains all information and settings git needs in order to run and version track your files. This also means that your git-tracked directory is self-contained, i.e. you can simply delete it and everything that has to do with git in connection to that directory will be gone.
* Notice that git told you: *nothing to commit (create/copy files and use "git add" to track*? Lets do that!  
Copy the following files and directories from the `reproducible_research_course/git_jupyter_docker` directory, into your `git_tutorial` directory:
  * `Dockerfile`
  * `Snakefile`
  * `config.yml`
  * `environment.yml`
  * `code/`  
* Once you have done that, run `git status` again. It will tell you that there are files in the directory that are not version tracked by git.

**Quick recap:** we have used two `git` commands this far:
* `git clone` - to clone a remote (e.g. online Bitbucket) repository locally. This is only run the first time one wants to download the repository locally.
* `git status` - this is a command you should use *a lot*. It will tell you, amongst other things, the status of your git clone in relation to the online remote repository.

### Committing

1. We will now commit the untracked files. A commit is essentially a set of changes to a set of files. Preferably, the changes making out a commit should be related to something, e.g. a specific bug fix or a new feature. Our first commit will be to add the copied files to the repository. Run (as suggested by `git status`):
```
git add Dockerfile Snakefile
```
* Run `git status` again! See that we have added Dockerfile and Snakefile to our upcoming commit (listed under *Changes to be committed*). This is called the staging area, and the files there are staged to be committed.
* We might as well commit all files in one go! So use `git add` on the remaining files as well:
```
git add config.yml environment.yml code/
```
Run `git status` and see that all files are in the staging area, and that no files are listed as untracked.
* We are now ready to commit! Run:
```
git commit -m "add initial files"
```
The `-m` option adds a commit message. This should be a short description of what the commit contains. If you forget to add `-m` and just run `git commit`, a terminal editor will open and prompt you to write a commit message. This can be confusing if you are note used to use a terminal editor, so try to remember the `-m` flag.
* Run `git status` (yep, again!). It should tell you *nothing to commit, working directory clean*.
* Now, let's edit a file. Open up `environment.yml` in your favorite editor, and change the version of bowtie2 to a different number, e.g. `bowtie2=2.1`.
* Run `git status`. It will tell you that there are modifications in one file (`environment.yml`) compared to the previous commit. This is nice! We don't have to keep track of what files we have edited, git will do that for us.
* Run `git diff environment.yml`. This will show you the changes made to the file. A `-` means a deleted line, a `+` means an added line. There are also shown a few lines before and after the changes, to put them in context.
* Let's edit another file! Open `config.yml`and change the line `genome_id: NCTC8325` to `genome_id: ST398`. Run `git status`. Run `git diff`. If we don't specify a file, it will show all changes made in any file, compared to the previous commit. Do you see your changes?
* Ok, we made our changes. Let's commit them! Run:
```
git add -A
```
This will add both our files to the staging area at the same time. Run `git status` and see that the changes in both `config.yml` and `environment.yml` are ready to be committed.
* But wait a minute! Shouldn't each commit optimally be a specified set of changes? Yes! So we want to make two commits, one for each change. Let's remove `environment.yml` from the staging area. `git status` tells us how to do this: *(use "git reset HEAD <file>..." to unstage)*. So run:
```
git reset HEAD environment.yml
```
Run `git status` again. See that now only `config.yml` is staged for being committed, whereas the changes in `environment.yml` are tracked by git, but not ready to be commited.
* Commit the changes in `config.yml`:
```
git commit -m "change genome to align to"
```
* Add and commit the changes in `environment.yml`:
```
git status
git add environment.yml
git status
git commit -m "change bowtie2 version"
git status
```
You don't have to run `git status` between each command, but it can be useful in the beginning while learning what each command does.
* To see a history of our changes so far, run:
```
git log
```
As you can see, each commit is a point in history. The more often you commit, and the more specific you keep your commits, the better (more fine-grained) history and version tracking you will have of your files.
* We can also try to delete a file:
```
rm Dockerfile
```
Run `git status`. As you can see, git tells us that the file is deleted, but the deletion is not committed. In the same way as we commit or edits to files, we need to commit a deletion of a file:
```
git rm Dockerfile
git status
git commit -m "remove Dockerfile"
git status
git log
```

**Quick recap:** we now added four important git commands to our repertoire:
* `git add` - add a file so that changes in that file can be committed.
* `git rm` - same as `git add` but sets a file to be deleted in the next commit.
* `git commit` - commits the changes we have staged (by using `git add` or `git rm`).
* `git log` - shows us the commit history.


### Pushing

So far we have just worked locally. The strength with git is that we can add a remote location to push our commits. In fact, we already have setup such a remote, since we created the repository at Bitbucket and cloned it locally. The idea is that you work and edit your files locally, and commit changes as you go along. At some points, preferably as often as possible, you push your changes to the remote, in this case Bitbucket. Your local copy and the remote copy are then in sync. In principle you can now safely delete your local copy and now that it is all backed up in the cloud, with the full commit history. This also enables collaboration. Several users can work on their local clones of a given repository and push changes to a common remote location. Let's try this out in practice!
1. Run `git remote -v`. This will show you what remote location is connected to your local git clone. The short name of the default remote is usually origin.
* Run `git branch`. This will show you the name of the current branch, by default this will be master.  
*We have not mentioned branches, and will not go into details about branches here, but they are a major feature of git. They allow you to have different "versions" of a repository. As an example, during software development it is common to have a release branch containing code that is working correctly, and a development branch containing code with new features and fixes but also potential bugs that have not been fixed yet. Once the development branch is fixed and working, it can be merged into the release branch. End-users will typically use the code in the release branch only. In the reproducible research setting, it is often enough to only use the master branch, but keep in mind that you have the possibility to use more branches if you see the needed.*
* Now we will push the latest commits to the master branch to our remote origin:
```
git push -u origin master
```
* Run `git status`. This should tell you that *Your branch is up-to-date with 'origin/master'.*.
* Go to your Bitbucket repository in your browser again and click on Source to the left. You should now see that the files you have locally appear here as well!
* Click on config.yml. You will see the contents of the file. Notice that it is the latest version, where we changed the genome_id.
* Click on Diff, in the upper left corner. You will see the changes made to this file compared to the previous commit.
* Click History. You will see an overview of the commits involving changes made to this file.
* Click Commits, to the left in the main menu. You will see an overview of all commits made. Click on a specific commit to see the changes introduced by that commit.
* Click on the commit that was the initial commit, where we added all the files. You can now click on View source in the top right corner. You will now see the files as they were when we first added them. Specifically you can see that the `Dockerfile` is back, even though we deleted it! Click on Source to the left again to return to the latest version.

**Quick recap:** we now learnt yet another important git command:
* `git push` -  to push local commits to the remote repository


### Conflicts

We will now learn how to manage conflicts. This is important to know, since it will probably happen sooner or later. The important thing is to not panic! :)

1. On the Bitbucket webpage, click on `environment.yml` and click Edit. We can now edit this file directly on the web. This is not mainly recommended, but we will do it to demonstrate a point.
* Let's pretend that using multiqc version 1.3 did not work. Change the multiqc version to 1.4:  
```
multiqc=1.4
```
Click Commit. Add the commit message: "update multiqc version to 1.4". Click Commit.
* Click Commits to the left to see the commit history, and your latest change at the top.
Now we have a change in the remote repository that is not yet in our local clone. This could happen for instance if a collaborator of yours committed a change and pushed it to Bitbucket.
* Go back to your local terminal. Run `git status`. Notice that git says: *Your branch is up-to-date with 'origin/master'.* This is of course not true, but our local git clone is not yet aware of the remote changes. We will get those changes soon.
* But first, we will edit `environment.yml` locally as well! (It may be the case that your collaborator thought it was good to use multiqc version 1.4, whereas you thought it would be better to use multiqc version 1.2, but neither of you communicated that to the other.) Use a text editor and change the multiqc line to:
```
multiqc=1.2
```
* Commit your change (use `git status` along the way if you want to check what is happening in the different steps):
```
git status
git add environment.yml
git status
git commit -m "downgraded multiqc to v1.2"
git status
```
* Now let's try to push this commit!
```
git push
```
Note that after the initial push you probably don't have to specify `-u origin master`, git will figure that out by itself. Read the error message. It should be fairly informative of what is going on. In essence it will not allow you to push since there are conflicting changes made to the remote repository.
* We will now dowload the changes made to the remote:
```
git fetch
```
* Now run `git status`. Unlike before, our local git clone now is aware of the latest changes pushed to the remote. It will tell you something along the lines: *Your branch and 'origin/master' have diverged,
and have 1 and 1 different commit each, respectively.*
* Now, since we ran `git fetch` our local git has up-to-date information about the status of the remote repository. We can therefore run
```
git diff origin/master
```
to see what the difference is between the current state of our local clone and the master branch on the remote origin, i.e. Bitbucket.
* Now let's try to integrate the remote changes with our local changes and get up to sync with the remote:
```
git pull
```
Note that you can skip the `git fetch`command if you want to and run `git pull` directly. The difference is that `fetch` will just update git with the latest information of the remote status, whereas `pull` will try to integrate and sync those changes to your local clone directly.
* As you have probably noticed, the `git pull` command resulted in a conflict. Git tells us about this and suggests that we should fix the conflicts and commit that.
As always, run `git status` to get an overview! You will see that you have, so called, unmerged paths and that the conflicting file is `environment.yml`, since both modified this file. To fix a conflict, open the affected file in a text editor. You will see that it now looks like this:
```yaml
  channels:
  - conda-forge
  - bioconda
  dependencies:
  - fastqc=0.11
  - sra-tools=2.8
  - bowtie2=2.1
  <<<<<<< HEAD
  - multiqc=1.2
  =======
  - multiqc=1.4
  >>>>>>> d9b35ef61d2fde56fcbd64aacb10a96098c67cbf
  - snakemake=4.3.0
  - samtools=1.6
  - htseq=0.9
  - graphviz=2.38.0
  - xorg-libxrender
  - xorg-libxpm
```
The part between `<<<<<<< HEAD` and `=======` is your local version, and the part between `=======` and `>>>>>>> d9b35ef61d2fde56fcbd64aacb10a96098c67cbf` is the one added to the remote and which caused the conflict when you tried to pull those changes to your local repository. The long sequence of numbers is the commit id (the first 7 are e.g. displayed on Bitbucket under Commits).
* It is now up to you to decide which version to keep, or to change it to a third alternative. Let's say that you are confident that it is better to run multiqc v1.2 rather than v1.4. Edit the file so that it looks like you want it to, i.e. remove the files added by git and delete the line with `multiqc=1.4`, the final file should look like this:
```yaml
  channels:
  - conda-forge
  - bioconda
  dependencies:
  - fastqc=0.11
  - sra-tools=2.8
  - bowtie2=2.1
  - multiqc=1.2
  - snakemake=4.3.0
  - samtools=1.6
  - htseq=0.9
  - graphviz=2.38.0
  - xorg-libxrender
  - xorg-libxpm
```
* Run `git status`, notice that it says *use "git add <file>..." to mark resolution*? Let's do that!
```
git add environment.yml
```
* Run `git status` again! It will now tell us: *All conflicts fixed but you are still merging. (use "git commit" to conclude merge)*. So, you probably guessed it, run:
```
git commit -m "merge and set multiqc to v1.2"
```
* Finally, push these changes to Bitbucket:
```
git push
```
* Go to Bitbucket in the browser and click Commits. You should be able to see a graph showing that the paths diverged (where one commit set the version to 1.4 and the other to 1.2) and that they are later merged, and the conflict fixed!

**Quick recap:** we now learnt how to sync or local clone with the remote one on Bitbucket, and how to fix potential conflicting commits. We added these commands to our repertoire:
* `git fetch` - downloads information from the remote repository.
* `git pull` - both fetches and integrates changes from the remote repository.

### Tagging

The last topic we will cover in this tutorial is tagging. Git allows us to tag commits. This is of particular convenience when thinking about reproducible research. We can tag commits that represent important points in history, like e.g. the version of the repository that was used for the submission, the version used during resubmission, and, most importantly, the version used for the final publication. The first two examples are mainly useful internally, but the latter is essential for other researchers being able to rerun your published analysis (assuming that you may continue to work on your code and push updates to the repository). Let's try this out!

1. Let's assume that the status of the repository as it is now is ready for a submission to a journal, e.g. it contains scripts that are used to generate the manuscript figures. Let's add a tag:
```
git tag "submission1"
```
* To push this tag to Bitbucket we use:
```
git push --tags
```
* Go to Bitbucket and check Commits. Can you see that the tag has been added?
* Let's assume we now got review comments and by fixing those had to update our code. Open `config.yml` and change the line `max_reads: 25000` to `max_reads: 50000`. Commit and push that change:
```
git add config.yml
git commit -m "increased number of reads"
git push
```
* Now let's say that the reviewers were happy and the manuscript was accepted for publication. Let's immediately add a tag:
```
git tag "publication"
git push --tags
```
* After the study was published you realized that you get nicer QC information if you upgrade multiqc. Open `environment.yml` and change `multiqc=1.2` to `multiqc=2.0`. Add, commit, and push this change:
```
git add environment.yml
git commit -m "upgrade to newer multiqc version"
git push
```
* Go to Bitbucket and click Downloads, and the Tags. Here users can download a compressed file containing the repository at the versions specified by the tags.
Alternatively, git users who want to reproduce your analysis with the code used for the publication can clone the Bitbucket repository and then run `git checkout publication`.
* You can try this in your local clone, run:
```
git checkout publication
```
Open `environment.yml` and note that the multiqc version now is back to 1.2! To go back to the latest version, run:
```
git checkout master
```
Again, open `environment.yml` and see that it now has version 2.0!

## Where to go next?

Now that you have learnt how to version track your code, you can continue to any other tutorial that will help you making your research even more reproducible.
