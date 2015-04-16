push:
	git commit -a
	git log --pretty --date=short > ChangeLog
	git commit --all --amend --no-edit
	git push
