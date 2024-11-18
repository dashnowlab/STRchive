import { truncate } from "lodash-es";
import { EChart } from "@kbox-labs/react-echarts";
import classes from "./RangeChart.module.css";

const RangeChart = ({
  title,
  xLabel,
  yLabels,
  values,
  counts,
  lowerBounds,
  upperBounds,
}) => {
  const upperMinusLower = lowerBounds.map(
    (_, index) => upperBounds[index] - lowerBounds[index],
  );

  return (
    <EChart
      height={50 + values.length * 40}
      className={classes.chart}
      grid={{ containLabel: true }}
      renderer="svg"
      textStyle={{ fontFamily: "inherit" }}
      title={{
        text: title,
        left: "center",
        textStyle: {
          fontSize: 20,
        },
      }}
      xAxis={{
        type: "value",
        name: xLabel,
        axisLabel: {
          fontSize: 20,
        },
        nameTextStyle: {
          fontSize: 20,
        },
        nameLocation: "middle",
        nameGap: 50,
      }}
      yAxis={{
        type: "category",
        data: yLabels,
        axisLabel: {
          fontSize: 20,
          padding: 20,
          formatter: (value) =>
            truncate(value, {
              length: window.innerWidth / 30,
            }),
        },
        nameTextStyle: {
          fontSize: 20,
        },
      }}
      series={[
        {
          data: lowerBounds,
          type: "bar",
          stack: "group",
          itemStyle: {
            opacity: 0,
          },
        },
        {
          data: upperMinusLower,
          type: "bar",
          barWidth: 10,
          stack: "group",
          itemStyle: {
            color: "var(--primary)",
            opacity: 1,
          },
          emphasis: {
            itemStyle: {
              color: "color-mix(in srgb, var(--primary), var(--white) 50%)",
            },
          },
        },
        {
          data: values,
          type: "scatter",
          itemStyle: {
            color: "var(--secondary)",
            opacity: 1,
          },
          emphasis: {
            itemStyle: {
              color: "color-mix(in srgb, var(--secondary), var(--white) 50%)",
            },
          },
          symbol: "diamond",
          symbolSize: 20,
        },
      ]}
      tooltip={{
        trigger: "axis",
        axisPointer: {
          type: "shadow",
        },
        formatter: ([lower, upperMinusLower, value]) =>
          `<div class="${classes.tooltip}">
            <p>
              <b>${value.name}</b>
            </p>
            <p>
              <span>Value:</span><br/>
              <b>${value.value.toFixed(2)}%</b>
            </p>
            <p>
              <span>95% Confidence Interval:</span><br/>
              <b>${lower.value.toFixed(2)}% &ndash; ${(lower.value + upperMinusLower.value).toFixed(2)}%</b>
            </p>
            <p>
              <span>Samples:</span><br/>
              <b>${counts[lower.dataIndex].toLocaleString()}</b>
            </p>
          </div>`,
      }}
    />
  );
};

export default RangeChart;
