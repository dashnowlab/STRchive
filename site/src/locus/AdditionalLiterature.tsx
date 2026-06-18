import type { Citation as CitationType } from "@/data/types";
import ShowMoreItems from "@/components/ShowMoreItems";
import Citation from "./Citation";

type Props = {
  additionalLiterature: CitationType[];
};

/** additional literature section */
export default function AdditionalLiterature({ additionalLiterature }: Props) {
  return (
    <ShowMoreItems
      items={additionalLiterature.map((citation, index) => (
        <Citation key={index} {...citation} />
      ))}
    />
  );
}
