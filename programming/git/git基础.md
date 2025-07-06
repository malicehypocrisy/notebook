```git
#查看配置
git config --list
#查看日志
ssh -v -i "D:/typora/Git/.ssh/id_rsa" -T git@github.com
#连接github
ssh -i "D:/typora/Git/.ssh/id_rsa" -T git@github.com
#初始化仓库
git init
#配置用户信息（只有配置一次）
git config --global user.name "你的GitHub用户名"
git config --global user.email "你的GitHub注册邮箱"
#关联本地仓库到远程
git remote add origin ssh
#删除已经存在URL
git remote remove origin
#提交文件
git add .          # 添加所有文件到暂存区
git commit -m "首次提交：添加全部笔记"  # 提交更改
#推到github
git push -u origin master # 首次推送需要加 `-u`
#日常推送
git add .
git commit -m "更新笔记：添加新内容"
git push


```

