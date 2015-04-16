push:
	git commit -a
	git log --pretty --date=short > ChangeLog
	git commit --amend --no-edit
	git push
