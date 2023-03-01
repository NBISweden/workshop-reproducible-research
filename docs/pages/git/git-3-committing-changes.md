We will now commit the untracked files. A commit is essentially a set of
changes to a set of files. Preferably, the changes making out a commit should
be related to something, *e.g.* a specific bug fix or a new feature.

* Our first commit will be to add the copied files to the repository. Run the
  following (as suggested by `git status`):

```bash
git add Dockerfile Snakefile
```

* Run `git status` again! See that we have added Dockerfile and Snakefile to
  our upcoming commit (listed under "*Changes to be committed*"). This is
  called the staging area, and the files there are staged to be committed.

* We might as well commit all files in one go! Use `git add` on the remaining
  files as well:

```bash
git add config.yml environment.yml
```

* Run `git status` and see that all files are in the staging area, and that no
  files are listed as untracked.

* We are now ready to commit! Run the following:

```bash
git commit -m "Add initial files"
```

The `-m` option adds a commit message. This should be a short description of
what the commit contains.

> **Good commit messages** <br>
> Writing informative and succinct commit messages can be tricky when you're
> just starting out. Here are some general guidelines that can help you write
> good commit messages from the start:
>
> * Separate subject from body with a blank line
> * Limit the subject line to 50 characters
> * Capitalize the subject line
> * Do not end the subject line with a period
> * Use the [imperative mood](https://en.wikipedia.org/wiki/Imperative_mood)
>   in the subject line
> * Wrap the body at 72 characters
> * Use the body to explain *what* and *why* vs. *how*
>
> In the command above we just added a short subject line ("Add initial
> files"). It is capitalized, less than 50 characters, does not end with
> a period, and uses imperative mood (Add!). It is possible to add
> a descriptive body text as well, as hinted by the points above. This is
> easiest done in a text editor. If you run `git commit` without the `-m`
> flag, Git will open the default terminal text editor (which can be
> configured with the `core.editor` variable) where you can write a longer
> commit message and body. If you want to read more about the motivation for
> these points, please see [this website](https://chris.beams.io/posts/git-commit/).

* Run `git status` again. It should tell you *"nothing to commit, working
  directory clean"*.

What have we done, so far? We had some files in our working directory that we
added to the Git staging area, which we subsequently committed to our Git
repository. A schematic overview of this process can be seen in the following
figure:

![](images/git_overview_local.png){ width=600px }

Let's repeat this process by editing a file!

* Open up `environment.yml` in your favorite editor, and change the version of
  bowtie2 to a different value, *e.g.* `bowtie2=2.2.4`.

* Run `git status`. It will tell you that there are modifications in one file
  (`environment.yml`) compared to the previous commit. This is nice! We don't
  have to keep track of which files we have edited, Git will do that for us.

* Run `git diff environment.yml`. This will show you the changes made to the
  file. A `-` means a deleted line, a `+` means an added line. There are also
  shown a few lines before and after the changes, to put them in context.

* Let's edit another file! Open `config.yml` and change the line `genome_id:
  NCTC8325` to `genome_id: ST398`. Run `git status`. Run `git diff`. If we
  don't specify a file, it will show all changes made in any file, compared to
  the previous commit. Do you see your changes?

* Ok, we made our changes. Let's commit them! Run:

```bash
git add config.yml environment.yml
```

This will add both our files to the staging area at the same time. Run `git
status` and see that the changes in both `config.yml` and `environment.yml` are
ready to be committed.

But wait a minute! Shouldn't each commit optimally be a conceptual unit of
change? Here we have one change to the genome ID used for an analysis and one
change where another software version is specified: these should probably be
separate. We thus want to make two commits, one for each change.

* Let's remove `environment.yml` from the staging area. `git status` tells us
  how to do this: *"(use "git reset HEAD <file>..." to unstage)"*. So run:

```bash
git reset HEAD environment.yml
```

!!! Note
    Maybe you didn't see the same message as indicated above? Is Git telling you
    to use a `git restore` instead? This is another one of Git's newer and
    experimental commands, which aims to remove some confusion about what
    commands do what (as many have multiple functions). While we have opted to
    stick with the old and stable commands until the new commands are no longer
    considered experimental, you are very welcome to use `git restore` instead
    of `git reset` to unstage the file above!

* Run `git status` again. See that now only `config.yml` is staged for being
  committed, whereas the changes in `environment.yml` are tracked by Git, but
  not ready to be committed.

* Commit the changes in `config.yml`:

```bash
git commit -m "Change to ST398 for alignment"
```

* Add and commit the changes in `environment.yml`:

```bash
git status
git add environment.yml
git status
git commit -m "Change bowtie2 version"
git status
```

You don't have to run `git status` between each command, but it can be useful
in the beginning while learning what each command does.

As you can see, each commit is a point in history. The more often you commit,
and the more specific you keep your commits, the better (more fine-grained)
history and version tracking you will have of your files.

* We can also try to delete a file:

```bash
rm Dockerfile
```

* Run `git status`. As you can see, Git tells us that the file is deleted, but
  that the deletion is not committed. In the same way as we commit edits to
  files, we need to commit a deletion of a file:

```bash
git add Dockerfile
git status
git commit -m "Remove Dockerfile"
git status
git log
```

Here we used `rm Dockerfile` to delete the file and `git add Dockerfile` to
stage the deletion. You can also use `git rm Dockerfile` to do both these
operations in one step.

* To see a history of our changes so far, run:

```bash
git log
```

!!! Tip
    Since Git keeps track of changes in text, *e.g.* code and text-based
    documentation, there are some files which you should *not* commit. Examples
    of such files are file formats that are not text-based, *e.g.* Microsoft
    Word/Excel files or PDFs - although one might sometimes want to track one of
    these files regardless, such as when you have a static PDF report you
    received from a sequencing platform that's never going to change. Other
    files you shouldn't track are vary large text files, *e.g.* those larger
    than 50 MB.

!!! Success "Quick recap"
    We added four important Git commands to our repertoire:
    
    * `git add` adds a file to the staging area
    * `git commit` commits the changes we have staged
    * `git rm` is shorthand for `rm <file>; git add <file>`
    * `git log` shows us the commit history
