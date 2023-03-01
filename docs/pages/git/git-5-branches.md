One of the most useful features of Git is called *branching*. Branching allows
you to diverge from the main line of work and edit or update your code and
files (*e.g.* to test out a new analysis or some experimental feature) without
affecting your main work. If the work you did in the branch turns out to be
useful you can merge that back into your `main` branch. On the other hand, if
the work didn't turn out as planned, you can simply delete the branch and
continue where you left off in your main line of work. Another use case for
branching is when you are working in a project with multiple people. Branching
can be a way of compartmentalizing your team's work on different parts of the
project and enables merging back into the `main` branch in a controlled
fashion; we will learn more about this in the section about working remotely.

* Let's start trying out branching! We can see the current branch by running:

```bash
git branch
```

This tells us that there is only the `main` branch at the moment.

!!! Note "Main and Master"
    If your branch is called `master` instead of `main` that's perfectly fine as
    well, but do check out the Git section of the [pre-course setup](../course-information/pre-course-setup.md)
    for more details about the choice of default branch names.

* Let's make a new branch:

```bash
git branch test_alignment
```

* Run `git branch` again to see the available branches. Do you note which one
  is selected as the active branch?

* Let's move to our newly created branch using the `checkout` command:

```bash
git checkout test_alignment
```

!!! Tip
    You can create and checkout a new branch in one line with `git checkout -b
    branch_name`.

Let's add some changes to our new branch! We'll use this to try out a different
set of parameters on the sequence alignment step of the case study project.

* Edit the `Snakefile` so that the shell command of the `align_to_genome` rule
  looks like this (add the `--very-sensitive-local` option):

```bash
bowtie2 --very-sensitive-local -x $indexBase -U {input.fastq} > {output} 2> {log}
```

* Add and commit the change!

* To get a visual view of your branches and commits you can use the command:

```bash
git log --graph --all --oneline
```

It is often useful to see what differences exist between branches.
You can use the `diff` command for this:

```bash
git diff main
```

This shows the difference between the active branch (`test_alignment`) and
`main` on a line-per-line basis. Do you see which lines have changed between
`test_alignment` and `main` branches?

!!! Tip
    We can also add the `--color-words` flag to `git diff`, which instead
    displays the difference on a word-per-word basis rather than line-per-line.

!!! Note
    Git is constantly evolving, along with some of its commands. While the
    `checkout` command is quite versatile (it's used for more than just switching
    branches), this versatility can sometimes be confusing. The Git team thus
    added a new command, `git switch`, that can be used instead. This command is
    still experimental, however, so we have opted to stick with `checkout` for
    the course - for now.

Now, let's assume that we have tested our code and the alignment analysis is run
successfully with our new parameters. We thus want to merge our work into the
`main` branch. It is good to start with checking the differences between
branches (as we just did) so that we know what we will merge.

* Checkout the branch you want to merge into, *i.e.* `main`:

```bash
git checkout main
```

* To merge, run the following code:

```bash
git merge test_alignment
```

Run `git log --graph --all --oneline` again to see how the merge commit brings
back the changes made in `test_alignment` to `main`.

!!! Tip
    If working on different features or parts of an analysis on different
    branches, and at the same time maintaining a working `main` branch for the
    stable code, it is convenient to periodically merge the changes made to
    `main` into relevant branches (*i.e.* the opposite to what we did above).
    That way, you keep your experimental branches up-to-date with the newest
    changes and make them easier to merge into `main` when time comes.

* If we do not want to do more work in `test_alignment` we can delete that
  branch:

```bash
git branch -d test_alignment
```

* Run `git log --graph --all --oneline` again. Note that the commits and
  the graph history are still there? A branch is simply a pointer to a
  specific commit, and that pointer has been removed.

!!! Tip
    There are many types of so-called "branching models", each with varying
    degrees of complexity depending on the developer's needs and the number of
    collaborators. While there certainly isn't a single branching model that
    can be considered to be the "best", it is very often most useful to keep it
    simple. An example of a simple and functional model is to have a `main`
    branch that is always working (*i.e.* can successfully run all your code
    and without known bugs) and develop new code on feature branches (one new
    feature per branch). Feature branches are short-lived, meaning that they
    are deleted once they are merged into `main`.

!!! Success "Quick recap"
    We have now learned how to divide our work into branches and how to manage
    them:

    - `git branch <branch>` creates a new branch.
    - `git checkout <branch>` moves the repository to the state in which the
    specified branch is currently in.
    - `git merge <branch>` merges the specified branch into the current one.
