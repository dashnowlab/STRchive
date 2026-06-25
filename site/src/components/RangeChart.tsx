import { EChart } from "@kbox-labs/react-echarts";
import { truncate } from "lodash-es";

type Props = {
  title: string;
  xLabel: string;
  yLabels: string[];
  values: number[];
  lowerBounds: number[];
  upperBounds: number[];
  tooltip?: (index: number) => string;
};

/** chart that shows values as points and ranges as bars */
export default function RangeChart({
  title,
  xLabel,
  yLabels,
  values,
  lowerBounds,
  upperBounds,
  tooltip = () => "",
}: Props) {
  return (
    <EChart
      height={50 + values.length * 40}
      className="w-full"
      grid={{ containLabel: true }}
      renderer="svg"
      textStyle={{ fontFamily: "inherit", color: "var(--color-black)" }}
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
            truncate(value, { length: window.innerWidth / 30 }),
        },
        nameTextStyle: {
          fontSize: 20,
        },
      }}
      series={[
        {
          /** https://echarts.apache.org/en/option.html#series-boxplot.data */
          data: values.map((value, index) => [
            lowerBounds[index],
            value,
            value,
            value,
            upperBounds[index],
          ]),
          type: "boxplot",
          itemStyle: {
            borderWidth: 3,
            borderColor: "var(--color-primary)",
          },
        },
        {
          data: values,
          type: "scatter",
          itemStyle: {
            color: "var(--color-black)",
            opacity: 1,
          },
          emphasis: {
            itemStyle: {
              color: "var(--color-black)",
            },
          },
          symbolSize: 12,
        },
      ]}
      tooltip={{
        trigger: "axis",
        // eslint-disable-next-line -- no types available
        formatter: (params: any) =>
          `<div class="flex flex-col items-center gap-4 max-w-full">
            ${tooltip(params?.[0]?.dataIndex)}
          </div>`,
      }}
    />
  );
}
