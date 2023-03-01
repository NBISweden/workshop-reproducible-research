It is not uncommon to run into conflicts when you are trying to merge separate
branches, and it's even more common when you're working in a collaborative
setting with remote repositories. It'll happen sooner or later, even if you're
only working locally, so it's important to know how to deal with them! We'll
now introduce a conflict on purpose, which we can then solve.

* Remember that we have two separate local copies of the same repository? Let's
  go into the first one, `git_tutorial`, and change the MultiQC version in the
  `environment.yml` file:

```yaml
multiqc=1.8
```

* Add, commit and push your change to the remote.

Now we have a change in our remote and *one* of our local copies, but not in the
other. This could happen if a collaborator of yours committed a change and
pushed it to GitHub. Let's create a conflict!

* Move into your other local repository, `git_remote_tutorial`, which doesn't
  have the new change. Run `git status`. Notice that Git says: "*Your branch is
  up-to-date with 'origin/main'.*". We know that this is not true, but this
  local clone is not yet aware of the remote changes.

* Let's change the `environment.yml` file in this local repository as well, but
  to version 1.6, instead! It may be the case that your collaborator thought it
  was good to use MultiQC version 1.8, whereas you thought it would be better to
  use MultiQC version 1.6, but neither of you communicated that to the other.

* Add and commit your change and try to push the commit, which should give you
  an error message that looks like this:

```no-highlight
 ! [rejected]        main -> main (fetch first)
error: failed to push some refs to 'https://github.com/user/git_tutorial.git'
hint: Updates were rejected because the remote contains work that you do
hint: not have locally. This is usually caused by another repository pushing
hint: to the same ref. You may want to first integrate the remote changes
hint: (e.g., 'git pull ...') before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.
```

This error message is thankfully quite informative in regards to what is going
on and what might be done about it. In essence it will not allow you to push
to the remote since there are conflicting changes made to it.

* Let's download the changes made to the remote, but without trying to merge
  them directly. This can be done using the following command:

```bash
git fetch
```

!!! Note
    The `fetch` command is very similar to `pull` in that it downloads remote
    changes that are not present locally, but differs in that it doesn't try to
    merge them locally; `pull` both downloads and merges (unless there's
    a conflict, in which case it will tell you so and raise an error like the
    one above). You can thus skip `fetch` and just do `pull` straight away, if
    you prefer.

* Now run `git status`. Unlike before, our local Git clone now is aware of the
  latest changes pushed to the remote. It will tell you something along the
  lines: "*Your branch and 'origin/main' have diverged, and have 1 and
  1 different commit each, respectively.*".

* We can now run the following to see what the difference is between the current
  state of our local clone and the `main` branch on the remote origin:

```bash
git diff origin/main
```

* Now let's try to integrate the remote changes with our local changes and get
  up to sync with the remote:

```bash
git merge
```

Unsurprisingly, the `git merge` command resulted in a conflict. Git tells us
about this and suggests that we should fix the conflicts and commit that.

* As always, run `git status` to get an overview: you will see that you have
  so-called unmerged paths and that the conflicting file is `environment.yml`,
  since both modified the same line in this file. To fix a conflict, open the
  affected file in a text editor. You will see that it now looks something like
  this:

```no-highlight
channels:
  - conda-forge
  - bioconda
  - main
  - r
dependencies:
  - python=3.9.12
  - fastqc=0.11.9
  - sra-tools=2.10.1
  - snakemake=7.3.8
<<<<<<< HEAD
  - multiqc=1.6
=======
  - multiqc=1.8
>>>>>>> refs/remotes/origin/main
  - bowtie2=2.4.5
  - tbb=2020.2
  - samtools=1.15.1
  - subread=2.0.1
  - bedtools=2.29.2
  - wget=1.20.3
  - graphviz=3.0.0
  - r-base=4.1.3
  - r-ggplot2=3.3.5
  - r-reshape2=1.4.4
  - r-stringi=1.7.6
  - r-pheatmap=1.0.12
  - r-rmarkdown=2.13
  - r-r.utils=2.11.0
  - bioconductor-rtracklayer=1.54.0
  - bioconductor-geoquery=2.62.0
  - xorg-libxrender
  - xorg-libxpm
```

The part between `<<<<<<< HEAD` and `=======` is your local version, and the
part between `=======` and `>>>>>>> refs/remotes/origin/main` is
the one added to the remote and which caused the conflict when you tried to merge
those changes to your local repository. It is now up to you to decide which
version to keep, or to change it to a third alternative.

* Let's say that you are confident that it is better to run MultiQC 1.6 rather
  than 1.8. Edit the file so that it looks like you want it to, *i.e.* remove
  the lines added by Git and delete the line with `multiqc=1.8`. The final file
  should look like this:

```no-highlight
channels:
  - conda-forge
  - bioconda
  - main
  - r
dependencies:
  - python=3.9.12
  - fastqc=0.11.9
  - sra-tools=2.10.1
  - snakemake=7.3.8
  - multiqc=1.6
  - bowtie2=2.4.5
  - tbb=2020.2
  - samtools=1.15.1
  - subread=2.0.1
  - bedtools=2.29.2
  - wget=1.20.3
  - graphviz=3.0.0
  - r-base=4.1.3
  - r-ggplot2=3.3.5
  - r-reshape2=1.4.4
  - r-stringi=1.7.6
  - r-pheatmap=1.0.12
  - r-rmarkdown=2.13
  - r-r.utils=2.11.0
  - bioconductor-rtracklayer=1.54.0
  - bioconductor-geoquery=2.62.0
  - xorg-libxrender
  - xorg-libxpm
```

* Run `git status` again. Notice that it says `use "git add <file>..." to mark
  resolution`? Let's do that!

```bash
git add environment.yml
```

* Run `git status` again! It will now tell us: `All conflicts fixed but you are
  still merging. (use "git commit" to conclude merge)`. So, you probably
  guessed it, run:

```bash
git commit -m "Merge and set multiqc to v1.6"
```

* Finally, push these changes to GitHub:

```bash
git push
```

* Go to GitHub in the browser and click the commit tracker again. You will see
  a list of commits including where MultiQC was first changed to version `1.7`
  from our previous work, then to `1.8`, `1.6` and, finally, followed by a merge
  where the version was set to `1.6`.

!!! Note
    While the example we've used here is from a collaborative setting, conflicts
    also arise when you are working alone. They usually happen when you have
    several feature branches that you want to merge into `main` and you've
    forgot to keep all branches up-to-date with each other.

!!! Success "Quick recap"
    We learned about how conflicting commits can happen and how to deal with
    them by inspecting the affected files and looking for the source of the
    conflict.
