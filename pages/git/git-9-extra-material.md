The following extra material contains some more advanced things you can do with
Git and the command line in general, which is not part of the main course
materials. All the essential skills of Git are covered by the previous
sections; the material here should be considered tips and tricks from people
who use Git every day. You thus don't need to use these things unless you want
to, and you can even skip this part of the lesson if you like!

If you are interested in learning more about Git in general, here are some
reading tips for you:

* [Git cheat-sheet](https://education.github.com/git-cheat-sheet-education.pdf)
* [A simple Git guide](http://rogerdudler.github.io/git-guide/)
* [Resources to learn Git]( https://try.github.io/levels/1/challenges/1)
* [Git reference manual](https://book.git-scm.com/docs)

## Forking

When you want to work on an Open Source project that is available on *e.g.*
GitHub, you usually don't have permission to directly push code to the project's
repository - this is so that the project's maintainers are the only ones that
can directly change anything in their codebase. How do you then contribute to
projects that don't allow you to push your code to their repository? Simple: use
*forking*!

Forking is when you make your own copy of a repository on your GitHub account,
which you will then have permissions to change as you see fit. You can then
create pull requests from your fork to the original repository, rather than
pushing code to a new branch and making a pull request from that. Working with
forks just adds an additional step to the whole workflow: instead of being
"clone; code and commit changes on a new branch; push branch to remote; pull
request from branch" it becomes "fork; clone; code and commit changes; push code
to fork; pull request from fork".

You might also want to do a fork of a project simply because you want to have
your own copy of it as well, without ever having the intention of changing it.
This is, of course, perfectly fine as well, but do keep in mind that developers
are usually quite happy to incorporate new changes from contributors if they are
reasonable and fulfil a purpose and add functionality to the project. It is
quite common that you have a use-case the maintainer didn't think of before, and
that you've helped the project grow by contributing your code!

## Amending commits

Once in a while you'll have just commited something to your Git repo and
immediately remembered that you forgot to add something small, or perhaps you
saw an error somewhere. While you can certainly just add that and make a new
commit, wouldn't it be nicer if you could just make the change as if it was
already a part of the first commit? Well, you can! Just make the change, stage
it and the commit together with the `--amend` flag, like so:

```bash
git add <file>
git commit --amend
```

This will add the staged changes to the previous commit as if they had always
been there. Be careful, though! This will actually rewrite history, meaning that
it only works if you only amended local changes. If you had already pushed the
first commit to a remote repository you would run into trouble: you will be able
to make the amend without issue, but you'll get an error when you try to push
your new changes, since the remote already contains the *first* version of the
commit and can't simply rewrite what it already has.

Amending changes is thus a good way to fix small mistakes you realise you made
just after commiting them, as long as you only amend local changes!

## Rebasing

The `git rebase` command is an alternative to `git merge` in that it solves the
same problem: getting changes in one branch into another branch. We've already
gone through merging extensively, so how is rebasing different? Let's look at a
common case: a `feature-branch` which we want to get into the `main` branch.

![](images/git-rebase-1.png){ width=300px }

Recall that a merge creates a *merge commit*, something akin to `Merge branch
'feature-branch' into main` or similar. This is a *new* commit that didn't exist
that brings the changes on `feature-branch` into `main`, but it contains no
actual work itself. This is both a good and a bad thing: good, because merging
is a safe, *non-destructive* operation (it doesn't alter history); bad, because
it can make the history itself looks quite messy. These are the commands used
and what the history will look like afterwards:

```bash
git checkout main
git merge feature-branch
```

![](images/git-rebase-2.png){ width=400px }

(The commit with the dashed border is the merge commit.)

Rebasing, on the other hand does *not* create merge commits. Indeed, what rebase
does is to "re-base" one branch on the other, *i.e.* pretend that new changes
were done on a different base than what actually happened (hence the name).
Getting our `feature-branch` onto `main` using rebase actually entails two
separate steps: first the rebase itself, followed by a fast-forward merge:

```bash
git checkout feature-branch
git rebase main
```

![](images/git-rebase-3.png){ width=500px }

This step rebases our `feature-branch` on top of `main`, meaning that we pretend
that the commits on `feature-branch` were done based on the latest commits on
`main` - you can also think of it as moving the entire `feature-branch` to the
tip of the `main` branch. The commits with the dashed borders here indicate
*brand new* commits; rebasing can't somehow move the commits to the new base,
rather it has to "replay" those commits as if they were done on the new base.

```bash
git checkout master
git merge feature-branch
```

![](images/git-rebase-4.png){ width=500px }

We've now got our `feature-branch` commits onto `main` with a single, linear
history without any merge commits! We did have to rewrite history, though, when
we did the rebase itself. As with amending (see [above](#amending-commits)),
this is fine if we're only working locally, but we'll quickly run into trouble
if we try to rebase things that have already been pushed. We can rebase *on top
of* remote things, of course, since we're not changing any remote history, only
the local history. Be careful when you rebase!

## Rebasing as cleanup

If the above section felt scary, don't worry! There's another highly useful
use-case for `git rebase` that doesn't risk destroying any history, namely local
cleanup!

Let's imagine you've worked on your local `feature-branch` for some time, and
you have a number of commits on it. Some are highly related to each other and
might actually be better suited as a single commit. You've also spotted a
spelling error in one commit message, and realised that you missed important
information in another. We can actually solve all of these issues with an
*interactive rebase*! If you have 4 commits on your branch you can type the
following:

```bash
git rebase -i HEAD~4
```

The `-i` flag means *interactive*, while `HEAD~4` means 4 commits back from
`HEAD`. This will open your default text editor and give you a selection looking
something like this:

```no-highlight
pick 0abf162 First feature commit
pick befc682 A minor change on the first commit
pick c9d1426 A commit with an uncomplete commit message
pick 2e0cb97 A commit with a spelling mitake

# Rebase 879ddcc..0abf162 onto 879ddcc (4 commands)
#
# Commands:
# p, pick <commit> = use commit
# r, reword <commit> = use commit, but edit the commit message
# e, edit <commit> = use commit, but stop for amending
# s, squash <commit> = use commit, but meld into previous commit

(... more instructions ...)
```

The commits are ordered with the most recent one at the bottom. The commented
instructions (all of which are not shown here) show you what alternatives you
have to work with; all you have to do is to change the `pick` keyword next to
the commit hashes to whatever keyword you need from the list, save and exit.

In order to solve the toy example here we might decide that the four keywords
should be `pick`, `squash`, `reword` and `reword`, from top to bottom. Once
that's done simply save and exit, and another instance of your default text
editor will open for you to complete the specified changes. In the case above
we'd get two separate new instances where we can change the commit message -
these work the same as any normal commit.

Interactive rebasing is thus well-suited for fixing and cleaning of local
changes you have yet to push anywhere, even if you don't use rebasing as an
alternative to merging! This can make your Git history both cleaner and more
concise, which is great when you're collaborating with others.

## The reflog

We have shown many ways to work with Git and its various commands, and it
occasionally happens that errors are introduced - especially when you're not
careful with using `git commit --amend` or `git rebase` on remote changes. This
is where the *reflog* comes in. Think of the reflog as Git's "safety net": it
stores almost every change you make to a Git repository (regardless of whether
you commit the change) in a chronological manner. The following is an example of
what the output of the `git reflog` command might show:

```no-highlight
58deba6 HEAD@{0}: merge: feature-branch: Fast-forward
8c80c88 HEAD@{1}: checkout: moving from feature-branch to main
555544a HEAD@{2}: commit: feature development 2
4c92630 HEAD@{3}: commit: feature development 1
8c80c88 HEAD@{4}: checkout: moving from main to feature-branch
```

It show the most recent change at the top, notified by `HEAD@{0}`. We thus have
a merging of `feature-branch` into `main`, a checkout into `main`, two commits
on `feature-branch` and a checkout into `feature-branch` - reading it backwards
we get a chronological log of what has happened.

The reflog is incredibly useful for when you've lost something you later realise
you want to access again, such as when you've just used `git reset`. The reflog
might look like this, for example:

```no-highlight
bc3641f HEAD@{0}: reset: moving to HEAD~2
caf9321 HEAD@{1}: commit: More work on the feature
1bc36af HEAD@{2}: commit: Work on a new feature
```

We see two commits related to some new feature and a reset to `HEAD~2` (two
commits back from `HEAD`). If we realise that we actually liked the work we just
threw away we can move around in the reflog in a similar manner we do normal
commits:

```bash
git checkout HEAD@{1}
```

This will put us back to the state we were in before we used `git reset`. We
here refer to the reflog using the `HEAD@{N}` notation, which differs from the
usual `HEAD~N` notation so that it is clear if it is the commit history or the
reflog that is intended. While the reflog is hopefully not something you'll have
to use often it's quite useful to know it exists, if only to be able to search
the internet for more details regarding a problem you've encountered!

## Decorating your prompt

When you are working on the command line interface (CLI), you will usually have
some small pieces of information relating to your current directory, the name
of the computer or host you're working on, and so forth. You've probably
already seen your prompt while working with Git throughout this lesson, but
here's an example of what one might look like:

```no-highlight
erikfmbp:~/teaching/workshop-reproducible-research erik.fasterius $
```

The above prompt contains the name of the computer, a colon, the current
working directory, the username and a dollar-sign; it is stored in the
variable `PS1`. You can type `echo $PS1` to see what variables your prompt
is made up of; the above example contains `\h:\W \u\$`, where `\h` is the
hostname, `\W` the working directory and `\u` the username.

Some people like to also show the current branch on their prompt, thus avoiding
having to type `git branch` continuously. There are several ways you might do
this, and we're only presenting one of them here: a bash function.

```bash
git_branch() {
     git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}
```

This function does a number of things:

1. Ejects the error message from Git if the current directory isn't a part of a
   Git repository into `/dev/null` (_i.e._ into nothing).
2. Find the current branch by searching for a line that starts with `*` (*i.e.*
   the current branch) using the command line program `sed`.
3. Put the current branch into parentheses with a space before it.

We can then build our new prompt by adding this function into it:

```bash
# The first part of the old prompt
PS1='\h:\W \u'

# Add the Git branch
PS1=$PS1'$(git_branch)'

# Add the last part of the old prompt
PS1=$PS1' \$'
```

Now you should see the current Git branch on your prompt! The only problem now
is that this only works for your current session: once you restart your CLI
you'll have to re-define your prompt again. This can be circumvented, though.
What you need to do is to add the code defining your prompt into your so-called
bash profile: `~/.bash_profile`. Every time you load a new CLI session this
file is read and any code inside it is executed. You might already have this
file, so make sure you don't overwrite it!

## Bash aliases for git

Some Git commands are used over and over again when working with git, such as
`git status`. Some people like to have aliases (*i.e.* shortcuts) for these
common commands. Here is a small list of such aliases that you may find useful
or, even better, might inspire you to create your own! Add them to your
`~/.bash_profile` as above, so that they're available across sessions.

```bash
# Basic git commands
alias gb='git branch'
alias ga='git add'
alias gd='git diff'
alias gcm='git commit'
alias gp='git push'
alias gu='git pull'
alias gm='git merge'
alias gco='git checkout'
alias gl='git log'

# Git status in short format
alias gst='git status -s'

# Add and commit all tracked and modified files
alias gca='git commit -a'

# Create and checkout a new branch
alias gcob='git checkout -b'

# Git log with one line per commit
alias glo='git log --oneline'
```
