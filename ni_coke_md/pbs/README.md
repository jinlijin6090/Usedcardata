# PBS 提交（PBS Pro，`select=` 语法）

## 用法

```bash
cd ni_coke_md
vi pbs/pbs.conf                      # 第一次先检查：队列名、每节点核数、walltime、software 标签
bash pbs/submit.sh --dry-run         # 只打印 qsub 命令和每个作业包含哪些算例，不提交
bash pbs/submit.sh                   # 提交全部任务
bash pbs/submit.sh task1_flux        # 只提交某个任务
python campaign/campaign.py status   # 看进度（也可以 qstat -u $USER）
```

`submit.sh` **可以反复运行**：已完成的算例和已经在队列里（排队或运行中）的算例会自动跳过，
只提交剩下的。所以任何时候出了问题（作业被杀、自动重交失败、节点故障），在登录节点重新运行
`bash pbs/submit.sh` 就能恢复。

## 它做了什么

1. **打包**：每个作业占一个节点，最多同时跑 `CORES_PER_NODE / THREADS_PER_CASE`（32/4 = 8）个算例。
   算例按剩余时长从长到短排序后分组，长算例和长算例分在一起，避免节点空等。
   不满 8 个算例的作业只申请需要的核数（`ncpus = 算例数 × 4`）。
   如果你们集群要求整节点独占，设 `FULL_NODE=1`。
2. **walltime 和自动续交**：`WALLTIME` 留空时，自动读取队列的 `resources_max.walltime`；
   读不到就用 24 小时。每个算例在 walltime 前 `MARGIN_H` 小时之后，遇到下一个断点就干净退出，
   作业随后把没跑完的算例**自动重新 qsub**（最多 `MAX_RESUB` 次）。
   出错的算例（退出码不是 0 或 3）不会重交，作业日志里会列出来。
3. **task2 自动提交**：task2 需要 task1 中 50 ps/C 的 3 个算例的 `final.data`。
   第一次提交时 task2 会显示"waiting"；task1 的作业结束时如果发现输入已经齐了，就自动提交 task2（`AUTO_TASK2=1`）。
4. **环境**：作业里不需要 `module load`，也不需要 conda，全部用包里自带的环境。
   脚本会清掉 `PYTHONPATH`，并把自带库放在 `LD_LIBRARY_PATH` 最前面，不会和 VASP 的 oneAPI 或 anaconda 环境冲突。

## 文件

| 文件 | 说明 |
|---|---|
| `pbs.conf` | 集群相关设置（唯一需要改的地方） |
| `submit.sh` | 提交 / 补交 |
| `job.pbs` | 作业脚本本体（由 submit.sh 带参数提交，不要直接 qsub） |
| `../runs/.jobs/` | 每个作业的算例列表 `*.list` 和 PBS 输出 |

## 默认设置下大概会提交什么

（`bash pbs/submit.sh --dry-run` 可以看到准确列表）

- task1：9 个算例 → 2 个作业。第一个作业：3 个 50 ps/C、3 个 20 ps/C、2 个 10 ps/C；第二个作业：1 个 10 ps/C（4 核）。
- task3：12 个算例 → 2 个作业（8 个 + 4 个）。
- task2：3 个算例 → 1 个作业（12 核），task1 完成后自动提交。

最长的算例（50 ps/C）大约需要 34 小时墙钟；如果 walltime 是 24 小时，会自动续交一次。

## 可能需要根据你们集群调整的地方

- `-l software=...`：你的 VASP 脚本里有 `-l software=VASP`。如果集群要求每个作业都声明软件，
  在 `PBS_SOFTWARE` 里填管理员给的值（比如 `LAMMPS` 或 `other`）。不要求就留空。
- 如果**计算节点上不能运行 qsub**，自动续交会失败，作业日志里会提示。
  这种情况下在登录节点定期运行 `bash pbs/submit.sh` 即可（可以放进 crontab）。
- 项目号或账户：写进 `EXTRA_QSUB`，例如 `EXTRA_QSUB="-P myproject"`。
