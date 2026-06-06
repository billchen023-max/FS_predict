# =========================
#  顶刊版式（Nature-like）Shiny App  |  Stable for older bslib
#  - 不用 font_system()/font_google()（避免版本/联网问题）
#  - 输入项加单位（更像临床工具/投稿）
#  - KPI 三卡 + 右侧对齐 Gauge（消除大空白）
# =========================

suppressPackageStartupMessages({
  library(shiny)
  library(xgboost)
  library(caret)
  library(pROC)
  library(ggplot2)
})



suppressPackageStartupMessages({
  library(shiny)
  library(caret)
  library(ggplot2)
  library(bslib)
  library(shinycssloaders)
})

pos_label <- "DN"

# ====== 载入模型与配置 ======
model <- readRDS("FINAL_reduced6_xgb_model.rds")
features <- readRDS("FINAL_reduced6_features.rds")
default_threshold <- readRDS("FINAL_reduced6_threshold_youden.rds")

# ====== 风险分层（可按临床改）======
risk_band <- function(p){
  if (p < 0.20) return("低风险")
  if (p < 0.50) return("中风险")
  return("高风险")
}
band_color <- function(band){
  switch(band,
         "低风险" = "#2E7D32",
         "中风险" = "#B26A00",
         "高风险" = "#B00020",
         "#37474F")
}

# ====== Nature-like CSS（克制高级）======
nature_css <- "
:root{
  --bg: #f6f7fb;
  --card: #ffffff;
  --text: #111827;
  --muted: #6b7280;
  --line: rgba(17,24,39,.10);
  --shadow: 0 10px 30px rgba(0,0,0,.06);
  --radius: 18px;
}
body{ background: var(--bg); color: var(--text); }
h1,h2,h3,h4{ letter-spacing: .2px; }
.small-muted{ color: var(--muted); font-size: 0.92rem; }
.cardlike{
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  box-shadow: var(--shadow);
}
.hrline{ height: 1px; background: var(--line); margin: 14px 0; }
.mono{ font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; }

.btn-primary{
  border-radius: 14px;
  font-weight: 650;
  padding: 10px 16px;
}
.btn-primary:hover{ transform: translateY(-1px); }
.btn-outline-secondary{ border-radius: 14px; padding: 10px 16px; }

.valuebox{
  padding: 16px 18px;
  border-radius: var(--radius);
  border: 1px solid var(--line);
  background: linear-gradient(180deg, rgba(255,255,255,1), rgba(255,255,255,.92));
  box-shadow: var(--shadow);
  min-height: 96px;
}
.valuebox .label{ color: var(--muted); font-size: .9rem; }
.valuebox .value{ font-size: 1.65rem; font-weight: 800; margin-top: 6px; }
.valuebox .sub{ color: var(--muted); margin-top: 6px; font-size: .92rem; }
.kpi-row{ margin-bottom: 8px; }
"

# ====== Gauge（半圆弧）======
gauge_plot <- function(p, th){
  # 0~1 映射到 180°~0°
  ang <- pi * (1 - p)
  th_ang <- pi * (1 - th)
  
  arc_df <- data.frame(
    t = seq(pi, 0, length.out = 220),
    x = cos(seq(pi, 0, length.out = 220)),
    y = sin(seq(pi, 0, length.out = 220))
  )
  
  needle <- data.frame(
    x = c(0, 0.92*cos(ang)),
    y = c(0, 0.92*sin(ang))
  )
  
  th_pt <- data.frame(
    x = 0.98*cos(th_ang),
    y = 0.98*sin(th_ang)
  )
  
  ggplot() +
    geom_path(data = arc_df, aes(x, y), linewidth = 10, alpha = 0.18) +
    geom_path(data = arc_df, aes(x, y), linewidth = 2, alpha = 0.35) +
    geom_segment(data = needle,
                 aes(x = x[1], y = y[1], xend = x[2], yend = y[2]),
                 linewidth = 1.7) +
    geom_point(aes(x = 0, y = 0), size = 4) +
    geom_point(data = th_pt, aes(x, y), size = 3, shape = 4, stroke = 1.3) +
    annotate("text", x = -1.0, y = -0.12, label = "0", size = 4, alpha = 0.7) +
    annotate("text", x =  1.0, y = -0.12, label = "1", size = 4, alpha = 0.7) +
    annotate("text", x =  0, y = -0.34,
             label = sprintf("P(DN)=%.3f   |   Threshold=%.3f", p, th),
             size = 4, alpha = 0.8) +
    coord_fixed(xlim = c(-1.15, 1.15), ylim = c(-0.45, 1.15), clip = "off") +
    theme_void() +
    ggtitle("Predicted probability gauge (DN)") +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# =========================
# UI
# =========================
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly"
  ),
  tags$head(tags$style(HTML(nature_css))),
  
  # 顶部标题栏
  div(class = "cardlike", style = "padding:18px 20px; margin: 16px 0 18px 0;",
      div(style="display:flex; justify-content:space-between; align-items:flex-end; gap:12px;",
          div(
            h3("纤维蛋白鞘（DN）风险评分工具", style="margin:0; font-weight:900;"),
            div(class="small-muted",
                "Reduced-6 XGBoost • 输出 P(DN | 6 inputs) • 阈值可调")
          ),
          div(class="small-muted mono",
              paste0("Model: ", "FINAL_reduced6_xgb_model.rds"))
      )
  ),
  
  # 主体两列：左 Input，右 Output
  layout_columns(
    col_widths = c(4, 8),
    
    # ===== 左侧输入卡 =====
    div(class = "cardlike", style="padding:18px 18px;",
        h4("Input", style="font-weight:900; margin-top:0;"),
        div(class="small-muted", "输入 6 项指标后点击“计算风险”。缺失值会提示。"),
        div(class="hrline"),
        
        # 这里的单位写法：你可以按你实际数据表改（例如 mm、cm、s、×10^9/L、fL）
        layout_columns(
          numericInput("Distance.to.RA", "Distance to RA (mm)", value = NA, min = -9999, max = 9999),
          numericInput("Distance.to.clavicle", "Distance to clavicle (mm)", value = NA, min = -9999, max = 9999),
          col_widths = c(6, 6)
        ),
        numericInput("Distance.to.bronchus", "Distance to bronchus (mm)", value = NA, min = -9999, max = 9999),
        
        layout_columns(
          numericInput("APTT", "APTT (s)", value = NA, min = 0, max = 500),
          numericInput("MON",  "Monocyte count (×10^9/L)", value = NA, min = 0, max = 100),
          col_widths = c(6, 6)
        ),
        numericInput("MPV", "MPV (fL)", value = NA, min = 0, max = 50),
        
        div(class="hrline"),
        
        sliderInput("threshold", "判定阈值（DN 阳性）",
                    min = 0, max = 1,
                    value = default_threshold,
                    step = 0.001),
        
        div(style="display:flex; gap:10px;",
            actionButton("calc", "计算风险", class = "btn-primary"),
            actionButton("reset", "清空输入", class = "btn btn-outline-secondary")
        ),
        
        div(style="margin-top:10px;",
            uiOutput("warn_box"))
    ),
    
    # ===== 右侧输出卡 =====
    div(class = "cardlike", style="padding:18px 18px;",
        h4("Output", style="font-weight:900; margin-top:0;"),
        div(class="small-muted", "概率、风险分层与二分类结果显示如下。"),
        div(class="hrline"),
        
        # KPI 三卡
        layout_columns(
          div(class="valuebox kpi-row",
              div(class="label", "DN probability"),
              div(class="value mono", textOutput("p_big", inline = TRUE)),
              div(class="sub", "P(DN | 6 inputs)")
          ),
          div(class="valuebox kpi-row",
              div(class="label", "Risk band"),
              div(class="value", uiOutput("band_big")),
              div(class="sub", "Low / Intermediate / High")
          ),
          div(class="valuebox kpi-row",
              div(class="label", "Binary decision"),
              div(class="value", uiOutput("pred_big")),
              div(class="sub", "Based on threshold")
          ),
          col_widths = c(4,4,4)
        ),
        
        # ===== 版式关键：左说明 + 右侧对齐 Gauge =====
        layout_columns(
          div(
            div(class="small-muted", style="font-weight:700; margin-bottom:6px;", "Interpretation"),
            tags$ul(
              class="small-muted",
              tags$li("Needle: predicted probability P(DN)."),
              tags$li("× marker: current threshold (adjustable)."),
              tags$li("Use as an adjunct to clinical judgement.")
            ),
            div(class="small-muted mono", style="margin-top:8px;",
                "Note: threshold default = Youden optimal point (test set).")
          ),
          
          div(style="display:flex; justify-content:flex-end;",
              div(style="width:560px;",
                  withSpinner(plotOutput("gauge_plot", height = "260px"), type = 6)
              )
          ),
          col_widths = c(5, 7)
        ),
        
        div(class="hrline"),
        
        # 可折叠说明（tags$details）
        tags$details(
          tags$summary("解释说明（可直接用于论文/工具说明）"),
          tags$ul(
            tags$li("阳性（DN）= 发生纤维蛋白鞘；阴性（Control）= 未发生。"),
            tags$li("输出为预测概率：P(DN | 6项输入)。"),
            tags$li("阈值用于二分类判定；默认阈值来自 test 集 Youden 最优点。"),
            tags$li("本页面仅做风险提示，最终以临床判断为准。")
          )
        )
    )
  )
)

# =========================
# Server
# =========================
server <- function(input, output, session) {
  
  # 清空输入
  observeEvent(input$reset, {
    updateNumericInput(session, "Distance.to.RA", value = NA)
    updateNumericInput(session, "Distance.to.clavicle", value = NA)
    updateNumericInput(session, "Distance.to.bronchus", value = NA)
    updateNumericInput(session, "APTT", value = NA)
    updateNumericInput(session, "MON", value = NA)
    updateNumericInput(session, "MPV", value = NA)
  })
  
  # 一行输入数据（注意：列名必须与模型训练时一致！这里保持原列名不变）
  make_one_row <- reactive({
    data.frame(
      `Distance.to.RA`       = as.numeric(input$Distance.to.RA),
      `Distance.to.clavicle` = as.numeric(input$Distance.to.clavicle),
      `Distance.to.bronchus` = as.numeric(input$Distance.to.bronchus),
      `APTT`                 = as.numeric(input$APTT),
      `MON`                  = as.numeric(input$MON),
      `MPV`                  = as.numeric(input$MPV),
      check.names = FALSE
    )
  })
  
  validate_inputs <- function(df){
    if (any(!is.finite(as.matrix(df)))) return(FALSE)
    TRUE
  }
  
  pred_res <- eventReactive(input$calc, {
    df <- make_one_row()
    
    if (!validate_inputs(df)) {
      return(list(ok = FALSE, msg = "❌ 有输入为空或非数字：请补全 6 项数值后再计算。"))
    }
    
    # caret 二分类概率列名用 pos_label
    p <- predict(model, newdata = df, type = "prob")[, pos_label]
    p <- as.numeric(p)
    
    th <- input$threshold
    pred_label <- ifelse(p >= th, "预测 DN（阳性）", "预测 Control（阴性）")
    
    list(
      ok = TRUE,
      p = p,
      th = th,
      pred = pred_label,
      band = risk_band(p)
    )
  })
  
  output$warn_box <- renderUI({
    r <- pred_res()
    if (is.null(r)) {
      div(class="small-muted", "提示：输入完成后点击“计算风险”。")
    } else if (!isTRUE(r$ok)) {
      div(class="alert alert-warning", style="border-radius:14px; margin:0;",
          r$msg)
    } else {
      div(class="alert alert-success", style="border-radius:14px; margin:0;",
          "✅ 计算完成。")
    }
  })
  
  output$p_big <- renderText({
    r <- pred_res()
    if (is.null(r) || !isTRUE(r$ok)) return("--")
    sprintf("%.3f", r$p)
  })
  
  output$band_big <- renderUI({
    r <- pred_res()
    if (is.null(r) || !isTRUE(r$ok)) return(HTML("<span style='color:#6b7280;'>--</span>"))
    col <- band_color(r$band)
    HTML(sprintf("<span style='color:%s; font-weight:900;'>%s</span>", col, r$band))
  })
  
  output$pred_big <- renderUI({
    r <- pred_res()
    if (is.null(r) || !isTRUE(r$ok)) return(HTML("<span style='color:#6b7280;'>--</span>"))
    HTML(sprintf("<span style='font-weight:900;'>%s</span>", r$pred))
  })
  
  output$gauge_plot <- renderPlot({
    r <- pred_res()
    if (is.null(r) || !isTRUE(r$ok)) return(NULL)
    gauge_plot(r$p, r$th)
  })
}

shinyApp(ui, server)
