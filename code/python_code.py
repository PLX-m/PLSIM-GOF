# -*- coding: utf-8 -*-
import random
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader, random_split
from sklearn.preprocessing import StandardScaler


# ------------------------------------------------------------
# 0) 固定随机种子
# ------------------------------------------------------------
def set_all_seeds(seed):
    
    if seed is None:
        return

    seed = int(seed)

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

    try:
        torch.use_deterministic_algorithms(True, warn_only=True)
    except Exception:
        pass

    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


# ------------------------------------------------------------
# 1) MLP 模型
# ------------------------------------------------------------
class MLP(nn.Module):
    def __init__(self, input_dim, output_dim,
                 hidden_dims=[64, 32],
                 activation=nn.ReLU,
                 dropout_rate=0.0):
        super().__init__()
        layers = []
        last = input_dim
        for h in hidden_dims:
            layers += [
                nn.Linear(last, h),
                nn.LayerNorm(h),
                activation(),
                nn.Dropout(dropout_rate)
            ]
            last = h
        layers.append(nn.Linear(last, output_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def _fit_Y_with_internal_val(
    Z_train, Y_train,
    lr=2e-4, epochs=5000,
    hidden=[32, 16],
    patience=200, tol=1e-5,
    batch_size=32, dropout_rate=0.0,
    weight_decay=5e-4,
    device=None,
    verbose=False,
    seed=None
):
    
    set_all_seeds(seed)

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    Z_arr = np.asarray(Z_train, dtype=np.float32)
    Y_arr = np.asarray(Y_train, dtype=np.float32)
    if Y_arr.ndim == 1:
        Y_arr = Y_arr.reshape(-1, 1)

    Z_t = torch.tensor(Z_arr, dtype=torch.float32).to(device)
    Y_t = torch.tensor(Y_arr, dtype=torch.float32).to(device)

    # 内部 80/20 切 train/val
    n_total = len(Z_t)
    n_tr = int(0.8 * n_total)
    n_val = n_total - n_tr

    if n_tr == 0 or n_val == 0:
        raise ValueError(f"Too few samples for inner train/val split: n_total={n_total}")
    else:
        dataset = TensorDataset(Z_t, Y_t)

        split_gen = torch.Generator()
        if seed is not None:
            split_gen.manual_seed(int(seed) + 11)

        train_ds, val_ds = random_split(
            dataset,
            [n_tr, n_val],
            generator=split_gen
        )

        if verbose:
            print(f"[Y-InnerSplit] n_total={n_total}, n_tr={n_tr}, n_val={n_val}")

    loader_gen = torch.Generator()
    if seed is not None:
        loader_gen.manual_seed(int(seed) + 12)

    train_loader = DataLoader(
        train_ds,
        batch_size=batch_size,
        shuffle=True,
        generator=loader_gen
    )

    val_loader = DataLoader(
        val_ds,
        batch_size=len(val_ds),
        shuffle=False
    )

    # 网络: Z_std -> 1 维 Y
    input_dim = Z_arr.shape[1]
    modelY = MLP(input_dim=input_dim,
                 output_dim=1,
                 hidden_dims=hidden,
                 dropout_rate=dropout_rate).to(device)

    optY = optim.Adam(modelY.parameters(), lr=lr, weight_decay=weight_decay)
    lossf = nn.MSELoss()

    best_loss = np.inf
    best_stateY = None
    best_ep = -1
    no_imp = 0

    train_hist, val_hist = [], []

    for ep in range(epochs):
        modelY.train()
        epoch_train_loss = 0.0
        nb = 0

        for Zb, Yb in train_loader:
            optY.zero_grad()
            predY = modelY(Zb).view(-1, 1)
            lossY = lossf(predY, Yb)
            lossY.backward()
            optY.step()
            epoch_train_loss += lossY.item()
            nb += 1

        epoch_train_loss /= max(nb, 1)

        # 验证
        modelY.eval()
        with torch.no_grad():
            val_loss = 0.0
            for Zb, Yb in val_loader:
                predY = modelY(Zb).view(-1, 1)
                val_loss = lossf(predY, Yb).item()

        train_hist.append(epoch_train_loss)
        val_hist.append(val_loss)

        if verbose and (ep % 50 == 0 or ep == epochs - 1):
            print(f"[Y Epoch {ep:4d}] train_loss={epoch_train_loss:.6f} | val_loss={val_loss:.6f}")

        # early stop 只看 Y 的 val_loss
        if val_loss + tol < best_loss:
            best_loss = val_loss
            best_stateY = {k: v.detach().clone() for k, v in modelY.state_dict().items()}
            best_ep = ep
            no_imp = 0
        else:
            no_imp += 1
            if no_imp >= patience:
                if verbose:
                    print(f"[Y Early Stop] best_val={best_loss:.6f} @ epoch={best_ep}, stop at epoch={ep}")
                break

    if best_stateY is not None:
        modelY.load_state_dict(best_stateY)
        if verbose:
            print(f"[Y Load Best] epoch={best_ep}, best_val={best_loss:.6f}")

    hist = {
        "train": train_hist,
        "val": val_hist,
        "best_epoch": best_ep,
        "best_val": best_loss
    }
    modelY._loss_history = hist

    return modelY


def _fit_X_with_internal_val(
    Z_train, X_train,
    lr=1e-4, epochs=5000,
    hidden=[32, 16],
    patience=100, tol=1e-4,
    batch_size=32, dropout_rate=0.0,
    weight_decay=1e-3,
    device=None,
    verbose=False,
    seed=None
):
    
    set_all_seeds(seed)

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    Z_arr = np.asarray(Z_train, dtype=np.float32)
    X_arr = np.asarray(X_train, dtype=np.float32)

    Z_t = torch.tensor(Z_arr, dtype=torch.float32).to(device)
    X_t = torch.tensor(X_arr, dtype=torch.float32).to(device)

    # 内部 80/20 切 train/val
    n_total = len(Z_t)
    n_tr = int(0.8 * n_total)
    n_val = n_total - n_tr

    if n_tr == 0 or n_val == 0:
        raise ValueError(f"Too few samples for inner train/val split: n_total={n_total}")
    else:
        dataset = TensorDataset(Z_t, X_t)

        split_gen = torch.Generator()
        if seed is not None:
            split_gen.manual_seed(int(seed) + 21)

        train_ds, val_ds = random_split(
            dataset,
            [n_tr, n_val],
            generator=split_gen
        )

        if verbose:
            print(f"[X-InnerSplit] n_total={n_total}, n_tr={n_tr}, n_val={n_val}")

    loader_gen = torch.Generator()
    if seed is not None:
        loader_gen.manual_seed(int(seed) + 22)

    train_loader = DataLoader(
        train_ds,
        batch_size=batch_size,
        shuffle=True,
        generator=loader_gen
    )

    val_loader = DataLoader(
        val_ds,
        batch_size=len(val_ds),
        shuffle=False
    )

    # 网络: Z_std -> p 维 X
    input_dim = Z_arr.shape[1]
    output_dim = X_arr.shape[1]

    modelX = MLP(input_dim=input_dim,
                 output_dim=output_dim,
                 hidden_dims=hidden,
                 dropout_rate=dropout_rate).to(device)

    optX = optim.Adam(modelX.parameters(), lr=lr, weight_decay=weight_decay)
    lossf = nn.MSELoss()

    best_loss = np.inf
    best_stateX = None
    best_ep = -1
    no_imp = 0

    train_hist, val_hist = [], []

    for ep in range(epochs):
        modelX.train()
        epoch_train_loss = 0.0
        nb = 0

        for Zb, Xb in train_loader:
            optX.zero_grad()
            predX = modelX(Zb)
            lossX = lossf(predX, Xb)
            lossX.backward()
            optX.step()
            epoch_train_loss += lossX.item()
            nb += 1

        epoch_train_loss /= max(nb, 1)

        # 验证
        modelX.eval()
        with torch.no_grad():
            val_loss = 0.0
            for Zb, Xb in val_loader:
                predX = modelX(Zb)
                val_loss = lossf(predX, Xb).item()

        train_hist.append(epoch_train_loss)
        val_hist.append(val_loss)

        if verbose and (ep % 50 == 0 or ep == epochs - 1):
            print(f"[X Epoch {ep:4d}] train_loss={epoch_train_loss:.6f} | val_loss={val_loss:.6f}")

        # early stop 只看 X 的 val_loss
        if val_loss + tol < best_loss:
            best_loss = val_loss
            best_stateX = {k: v.detach().clone() for k, v in modelX.state_dict().items()}
            best_ep = ep
            no_imp = 0
        else:
            no_imp += 1
            if no_imp >= patience:
                if verbose:
                    print(f"[X Early Stop] best_val={best_loss:.6f} @ epoch={best_ep}, stop at epoch={ep}")
                break

    if best_stateX is not None:
        modelX.load_state_dict(best_stateX)
        if verbose:
            print(f"[X Load Best] epoch={best_ep}, best_val={best_loss:.6f}")

    hist = {
        "train": train_hist,
        "val": val_hist,
        "best_epoch": best_ep,
        "best_val": best_loss
    }
    modelX._loss_history = hist

    return modelX


def estimate_conditional_expectations_cv_with_folds(
    Z, Y, X, fold_id,
    verbose=True,
    seed=None
):
    
    set_all_seeds(seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[ECE-CV-FOLDS] device={device}")

    fold_id = np.asarray(fold_id).reshape(-1)
    labels = np.unique(fold_id)

    n, q = Z.shape
    p = X.shape[1]
    E_YZ = np.zeros(n, dtype=np.float64)
    E_XZ = np.zeros((n, p), dtype=np.float64)

    for k in labels:
        te = (fold_id == k)
        tr = ~te

        Z_tr, Y_tr, X_tr = Z[tr], Y[tr], X[tr]
        Z_te = Z[te]

        if verbose:
            print(f"\n[Fold {k}] train={tr.sum()} | test={te.sum()} | "
                  f"input_dim(Z)={q} | out_dim(Y)=1 | out_dim(X)={p}")

        # ==============================
        # 1) 只在 train fold 上对 Z 拟合 scaler
        # ==============================
        scZ = StandardScaler().fit(Z_tr)
        Z_tr_std = scZ.transform(Z_tr)
        Z_te_std = scZ.transform(Z_te)

        # Y_tr 和 X_tr 保持原始尺度，不进行标准化
        Y_tr = Y_tr.reshape(-1, 1) if len(Y_tr.shape) == 1 else Y_tr
        X_tr = X_tr.reshape(-1, p) if len(X_tr.shape) == 1 else X_tr

        # ==============================
        # 2) 训练 Y 的网络
        # ==============================
        modelY = _fit_Y_with_internal_val(
            Z_train=Z_tr_std,
            Y_train=Y_tr,
            verbose=False,
            seed=None if seed is None else int(seed) + 1000 + int(k)
        )

        # ==============================
        # 3) 训练 X 的网络
        # ==============================
        modelX = _fit_X_with_internal_val(
            Z_train=Z_tr_std,
            X_train=X_tr,
            verbose=False,
            seed=None if seed is None else int(seed) + 2000 + int(k)
        )

        # ==============================
        # 4) 在 test fold 上做预测
        # ==============================
        Z_te_t = torch.tensor(Z_te_std, dtype=torch.float32).to(device)

        modelY.eval()
        modelX.eval()
        with torch.no_grad():
            EY_te = modelY(Z_te_t).cpu().numpy().reshape(-1)
            EX_te = modelX(Z_te_t).cpu().numpy()

        E_YZ[te] = EY_te
        E_XZ[te, :] = EX_te

        if verbose:
            histY = getattr(modelY, "_loss_history", None)
            histX = getattr(modelX, "_loss_history", None)
            if histY is not None:
                print(f"[Fold {k}] Y best_val={histY['best_val']:.6f} "
                      f"@ epoch={histY['best_epoch']}")
            if histX is not None:
                print(f"[Fold {k}] X best_val={histX['best_val']:.6f} "
                      f"@ epoch={histX['best_epoch']}")

    return [E_YZ, E_XZ]



def fit_conditional_expectation_models(
    Z_train,
    Y_train,
    X_train,
    verbose=False,
    seed=None
):
    
    set_all_seeds(seed)

    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )

    # =============================================================
    # 1) standardize Z
    # =============================================================

    scZ = StandardScaler().fit(Z_train)
    Z_tr_std = scZ.transform(Z_train)

    # =============================================================
    # 2) fit Y model
    # =============================================================

    modelY = _fit_Y_with_internal_val(
        Z_train=Z_tr_std,
        Y_train=Y_train,
        verbose=verbose,
        seed=None if seed is None else int(seed) + 100
    )

    # =============================================================
    # 3) fit X model
    # =============================================================

    modelX = _fit_X_with_internal_val(
        Z_train=Z_tr_std,
        X_train=X_train,
        verbose=verbose,
        seed=None if seed is None else int(seed) + 200
    )

    # =============================================================
    # 4) return fitted objects
    # =============================================================

    out = {
        "scalerZ": scZ,
        "modelY": modelY,
        "modelX": modelX,
        "device": device
    }

    return out



def predict_conditional_expectations(
    fit_obj,
    Z_new
):

    scZ = fit_obj["scalerZ"]
    modelY = fit_obj["modelY"]
    modelX = fit_obj["modelX"]
    device = fit_obj["device"]

    # ============================================================
    # standardize new Z
    # ============================================================

    Z_std = scZ.transform(Z_new)

    Z_t = torch.tensor(Z_std, dtype=torch.float32).to(device)

    # ============================================================
    # prediction
    # ============================================================

    modelY.eval()
    modelX.eval()

    with torch.no_grad():
        EY_hat = modelY(Z_t).cpu().numpy().reshape(-1)
        EX_hat = modelX(Z_t).cpu().numpy()

    return EY_hat, EX_hat